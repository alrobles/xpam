# test-pam_polygons.R -- Phase D: PolygonOccSource (streaming rasteriser).
# Shared fixtures (rect_cells/rect_wkt/read_xpoc/make_world) live in
# helper-polygons.R.

test_that("polygon source rasterises to the expected (species, cell) pairs", {
  skip_without_gdal()
  skip_on_cran()
  w <- make_world()
  on.exit(unlink(c(w$fenv, w$fv)), add = TRUE)

  lv <- c("sp1", "sp2", "sp3", "sp4", "sp5")
  expect_warning(
    eng <- pam_env_engine(w$fenv, occ_polygons = w$fv,
                          poly_field = "sciname",
                          species = lv, mins = c(-3, -3), maxs = c(3, 3),
                          nbs = c(4, 5), tile_rows = 7),
    "skipped")
  tru <- xpam_truth(terra::values(w$r, mat = TRUE), w$occ, lv,
                    c(-3, -3), c(3, 3), c(4, 5))
  expect_equal(eng$richness, tru$richness)
  expect_equal(eng$footprint, tru$footprint)
  expect_equal(eng$nnz, tru$nnz)
  expect_equal(eng$availability, tru$availability)
  expect_equal(eng$species, lv)
  expect_equal(eng$n_skipped, 2)   # "unmatched" level + the NA-field feature
})

test_that("all_touched = FALSE drops the sub-pixel triangle", {
  skip_without_gdal()
  skip_on_cran()
  nr <- 20; nc <- 30
  r <- xpam_sim_env(nr, nc, 1L, seed = 2L)
  fenv <- tempfile(fileext = ".tif")
  fv <- tempfile(fileext = ".gpkg")
  terra::writeRaster(r, fenv, overwrite = TRUE)
  on.exit(unlink(c(fenv, fv)), add = TRUE)
  # Triangle wholly inside pixel (col 5, row 14): touches it, centres none.
  tri <- "POLYGON ((5.1 5.1, 5.4 5.1, 5.1 5.4, 5.1 5.1))"
  v <- terra::vect(data.frame(sciname = "sp1", wkt = tri), geom = "wkt")
  terra::writeVector(v, fv, overwrite = TRUE)

  at <- pam_env_engine(fenv, occ_polygons = fv, poly_field = "sciname",
                       species = "sp1", mins = -3, maxs = 3, nbs = 4,
                       all_touched = TRUE)
  ct <- pam_env_engine(fenv, occ_polygons = fv, poly_field = "sciname",
                       species = "sp1", mins = -3, maxs = 3, nbs = 4,
                       all_touched = FALSE)
  expect_equal(at$nnz, 1)                    # one cell touched
  expect_equal(ct$nnz, 0)                    # no pixel centre inside
  expect_equal(at$footprint, 1)
  expect_equal(ct$footprint, 0)
})

test_that("occ_dump writes a sorted XPOC reusable via occ_file", {
  skip_without_gdal()
  skip_on_cran()
  w <- make_world()
  fx <- tempfile(fileext = ".xpoc")
  on.exit(unlink(c(w$fenv, w$fv, fx)), add = TRUE)

  expect_warning(
    eng <- pam_env_engine(w$fenv, occ_polygons = w$fv,
                          poly_field = "sciname",
                          occ_dump = fx, mins = c(-3, -3), maxs = c(3, 3),
                          nbs = c(4, 5), poly_chunk_rows = 3),
    "skipped")
  expect_true(file.exists(fx))
  expect_equal(eng$n_skipped, 1)             # only the NA-field feature
  xp <- read_xpoc(fx)
  expect_equal(xp$n_species, 6)              # derived levels, "unmatched" in
  expect_setequal(eng$species,
                  c("sp1", "sp2", "sp3", "sp5", "sp4", "unmatched"))
  # the dumped table is exactly the expected pair set
  got <- paste(xp$sp, xp$cell)
  occ_codes <- match(w$occ$species, eng$species) - 1L
  want <- paste(occ_codes, w$occ$cell - 1L)
  # dumped table also covers the "unmatched" species' block
  want <- c(want, paste(match("unmatched", eng$species) - 1L,
                        rect_cells(0, 2, 0, 2, w$nc)))
  expect_setequal(got, want)
  expect_equal(length(got), length(unique(got)))

  # same marginals when streamed back through occ_file
  eng2 <- pam_env_engine(w$fenv, occ_file = fx,
                         mins = c(-3, -3), maxs = c(3, 3), nbs = c(4, 5))
  expect_equal(eng2$richness, eng$richness)
  expect_equal(eng2$footprint, eng$footprint)
  expect_equal(eng2$nnz, eng$nnz)
})

test_that("polygon SRS differing from the grid is reprojected", {
  skip_without_gdal()
  skip_on_cran()
  nr <- 18; nc <- 36
  r <- terra::rast(nrows = nr, ncols = nc, nlyrs = 1,
                   extent = c(-180, 180, -90, 90), crs = "EPSG:4326",
                   vals = runif(nr * nc, -2, 2))
  fenv <- tempfile(fileext = ".tif")
  terra::writeRaster(r, fenv, overwrite = TRUE)
  fv <- tempfile(fileext = ".gpkg")
  on.exit(unlink(c(fenv, fv)), add = TRUE)

  # lon/lat rect covering x (-9.5, 0.5) x y (10.5, 19.5); project to 3857
  wkt <- "POLYGON ((-9.5 10.5, 0.5 10.5, 0.5 19.5, -9.5 19.5, -9.5 10.5))"
  v4326 <- terra::vect(data.frame(sciname = "sp1", wkt = wkt),
                       geom = "wkt", crs = "EPSG:4326")
  v3857 <- terra::project(v4326, "EPSG:3857")
  terra::writeVector(v3857, fv, overwrite = TRUE)

  eng <- pam_env_engine(fenv, occ_polygons = fv, poly_field = "sciname",
                        species = "sp1", mins = -3, maxs = 3, nbs = 4)
  # ALL_TOUCHED: the rect's edges at x=-9.5/0.5 fall inside cols 17 and 18
  # ([-10,0) and [0,10)); y 10.5..19.5 lies wholly in row 7 ([10,20)).
  cells <- 7 * nc + c(17, 18)
  tru <- xpam_truth(cbind(terra::values(r)[, 1]),
                    data.frame(cell = cells + 1L, species = "sp1"),
                    "sp1", -3, 3, 4)
  expect_equal(eng$nnz, tru$nnz)
  expect_equal(eng$footprint, tru$footprint)
})

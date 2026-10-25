# test-pam_poly_multi.R -- generalised polygon inputs (issue: IUCN-style
# per-group datasets): vectors of paths, directories, SpatVector objects,
# and the k-way merged occurrence stream (MergeOccSource).
#
# rect_cells/rect_wkt/make_world/read_xpoc live in helper-polygons.R.

.make_two_files <- function(nr, nc) {
  fenv <- tempfile(fileext = ".tif")
  terra::writeRaster(xpam_sim_env(nr, nc, 2L, seed = 1L), fenv,
                     overwrite = TRUE)
  # File 1: sp1 block + half of sp2.  File 2: sp3 + the OTHER half of sp2
  # plus a duplicate of sp1's exact block (cross-dataset dedup check).
  mk <- function(df) {
    f <- tempfile(fileext = ".gpkg")
    terra::writeVector(terra::vect(df, geom = "wkt"), f, overwrite = TRUE)
    f
  }
  f1 <- mk(data.frame(
    sciname = c("sp1", "sp2"),
    wkt = c(rect_wkt(2, 6, 3, 8, nr),       # sp1: cols 2:6 rows 3:8
            rect_wkt(10, 12, 12, 18, nr)))) # sp2 left half
  f2 <- mk(data.frame(
    sciname = c("sp3", "sp2", "sp1"),
    wkt = c(rect_wkt(20, 25, 1, 4, nr),     # sp3
            rect_wkt(12, 14, 12, 18, nr),   # sp2 right half
            rect_wkt(2, 6, 3, 8, nr))))     # sp1 duplicated
  occ <- rbind(
    data.frame(cell = c(rect_cells(2, 6, 3, 8, nc)) + 1L, species = "sp1"),
    data.frame(cell = c(rect_cells(10, 14, 12, 18, nc)) + 1L,
               species = "sp2"),
    data.frame(cell = c(rect_cells(20, 25, 1, 4, nc)) + 1L,
               species = "sp3"))
  list(fenv = fenv, files = c(f1, f2), occ = occ, nr = nr, nc = nc)
}

test_that("two polygon datasets merge into one PAM with dedup", {
  skip_without_gdal()
  skip_on_cran()
  w <- .make_two_files(20, 30)
  on.exit(unlink(c(w$fenv, w$files)), add = TRUE)

  lv <- c("sp1", "sp2", "sp3")
  eng <- pam_env_engine(w$fenv, occ_polygons = w$files,
                        poly_field = "sciname", species = lv,
                        mins = c(-3, -3), maxs = c(3, 3), nbs = c(4, 5),
                        tile_rows = 7)
  tru <- xpam_truth(terra::values(terra::rast(w$fenv), mat = TRUE),
                    w$occ, lv, c(-3, -3), c(3, 3), c(4, 5))
  expect_equal(eng$richness, tru$richness)
  expect_equal(eng$footprint, tru$footprint)
  expect_equal(eng$nnz, tru$nnz)            # sp1 duplicated: counted once
  expect_equal(eng$species, lv)
})

test_that("a directory of vector files is scanned and merged", {
  skip_without_gdal()
  skip_on_cran()
  w <- .make_two_files(20, 30)
  d <- tempfile(); dir.create(d)
  ok <- file.copy(w$files, file.path(d, basename(w$files)))
  on.exit(unlink(c(w$fenv, w$files, d), recursive = TRUE), add = TRUE)

  lv <- c("sp1", "sp2", "sp3")
  eng <- pam_env_engine(w$fenv, occ_polygons = d,
                        poly_field = "sciname", species = lv,
                        mins = c(-3, -3), maxs = c(3, 3), nbs = c(4, 5))
  tru <- xpam_truth(terra::values(terra::rast(w$fenv), mat = TRUE),
                    w$occ, lv, c(-3, -3), c(3, 3), c(4, 5))
  expect_equal(eng$nnz, tru$nnz)
  expect_equal(eng$richness, tru$richness)
})

test_that("an in-memory SpatVector works like a file path", {
  skip_without_gdal()
  skip_on_cran()
  nr <- 20; nc <- 30
  fenv <- tempfile(fileext = ".tif")
  terra::writeRaster(xpam_sim_env(nr, nc, 2L, seed = 1L), fenv,
                     overwrite = TRUE)
  on.exit(unlink(fenv), add = TRUE)

  v <- terra::vect(data.frame(
    sciname = c("sp1", "sp2"),
    wkt = c(rect_wkt(2, 6, 3, 8, nr), rect_wkt(10, 14, 12, 18, nr))),
    geom = "wkt")
  lv <- c("sp1", "sp2")
  eng <- pam_env_engine(fenv, occ_polygons = v, poly_field = "sciname",
                        species = lv,
                        mins = c(-3, -3), maxs = c(3, 3), nbs = c(4, 5))
  occ <- rbind(
    data.frame(cell = rect_cells(2, 6, 3, 8, nc) + 1L, species = "sp1"),
    data.frame(cell = rect_cells(10, 14, 12, 18, nc) + 1L, species = "sp2"))
  tru <- xpam_truth(terra::values(terra::rast(fenv), mat = TRUE),
                    occ, lv, c(-3, -3), c(3, 3), c(4, 5))
  expect_equal(eng$nnz, tru$nnz)
  expect_equal(eng$richness, tru$richness)
  expect_equal(eng$footprint, tru$footprint)
})

test_that("multiple datasets without explicit species levels error", {
  skip_without_gdal()
  skip_on_cran()
  w <- .make_two_files(20, 30)
  on.exit(unlink(c(w$fenv, w$files)), add = TRUE)

  expect_error(
    pam_env_engine(w$fenv, occ_polygons = w$files,
                   poly_field = "sciname",
                   mins = c(-3, -3), maxs = c(3, 3), nbs = c(4, 5)),
    "species")
})

test_that("occ_dump over merged datasets writes the unified XPOC", {
  skip_without_gdal()
  skip_on_cran()
  w <- .make_two_files(20, 30)
  fx <- tempfile(fileext = ".xpoc")
  on.exit(unlink(c(w$fenv, w$files, fx)), add = TRUE)

  lv <- c("sp1", "sp2", "sp3")
  eng <- pam_env_engine(w$fenv, occ_polygons = w$files,
                        poly_field = "sciname", species = lv,
                        occ_dump = fx,
                        mins = c(-3, -3), maxs = c(3, 3), nbs = c(4, 5),
                        poly_chunk_rows = 3)
  xp <- read_xpoc(fx)
  expect_equal(xp$n_species, 3)
  # dump = exactly the merged pair stream: sorted, deduped
  want <- paste(match(w$occ$species, lv) - 1L, w$occ$cell - 1L)
  got <- paste(xp$sp, xp$cell)
  expect_setequal(got, want)
  expect_equal(length(got), length(unique(got)))

  # and it reads back identically
  eng2 <- pam_env_engine(w$fenv, occ_file = fx,
                         mins = c(-3, -3), maxs = c(3, 3), nbs = c(4, 5))
  expect_equal(eng2$nnz, eng$nnz)
  expect_equal(eng2$richness, eng$richness)
})

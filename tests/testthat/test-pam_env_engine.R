test_that("GDAL engine agrees with the in-memory pam_env path", {
  skip_without_gdal()
  skip_on_cran()

  nr <- 40; nc <- 30
  set.seed(11)
  r1 <- terra::rast(nrows = nr, ncols = nc, vals = runif(nr * nc, -3, 3))
  r2 <- terra::rast(nrows = nr, ncols = nc, vals = rnorm(nr * nc))
  f1 <- tempfile(fileext = ".tif"); f2 <- tempfile(fileext = ".tif")
  terra::writeRaster(r1, f1, overwrite = TRUE)
  terra::writeRaster(r2, f2, overwrite = TRUE)
  on.exit(unlink(c(f1, f2)), add = TRUE)

  occ <- data.frame(
    cell = sample.int(nr * nc, 400, replace = TRUE),
    species = sample(paste0("sp", 1:12), 400, replace = TRUE))

  eng <- pam_env_engine(c(f1, f2), occ, nbs = c(4, 5),
                        mins = c(-4, -5), maxs = c(4, 5),
                        tile_rows = 7, threads = 2)

  env <- data.frame(cell = seq_len(nr * nc),
                    x = terra::values(r1)[, 1], y = terra::values(r2)[, 1])
  # reference: same 4x5 grid via the in-memory kernel
  bins <- xpam:::c_bin_assign_d(cbind(env$x, env$y), c(-4, -5), c(4, 5),
                               c(4L, 5L))
  idx <- match(occ$cell, env$cell)
  bp <- bins[idx]
  sp <- match(occ$species, unique(occ$species)) - 1L
  refk <- xpam:::c_env_pam_bitmap(sp, bp, length(unique(occ$species)), 20L)

  expect_equal(eng$kernel, "bitmap")
  expect_equal(eng$nnz, refk$nnz)
  expect_equal(eng$richness, refk$richness)
  expect_equal(eng$footprint, refk$footprint)
  expect_equal(eng$ncols, nc); expect_equal(eng$nrows, nr)
  expect_gt(eng$n_strips, 1)                       # strip loop exercised
  # availability = cells per bin over ALL finite env cells
  av <- xpam:::c_availability(bins, 20L)
  expect_equal(eng$availability, as.numeric(av))
  expect_equal(eng$species, unique(occ$species))
})

test_that("engine falls back to the sort kernel past the budget", {
  skip_without_gdal()
  skip_on_cran()

  r <- terra::rast(nrows = 20, ncols = 20, vals = runif(400))
  f <- tempfile(fileext = ".tif")
  terra::writeRaster(r, f, overwrite = TRUE)
  on.exit(unlink(f), add = TRUE)

  occ <- data.frame(cell = 1:20, species = paste0("s", 1:20))
  res <- pam_env_engine(f, occ, nbs = 4L, mins = 0, maxs = 1,
                        max_bitset_mb = 0)
  expect_equal(res$kernel, "sort")
  big <- pam_env_engine(f, occ, nbs = 4L, mins = 0, maxs = 1)
  expect_equal(big$kernel, "bitmap")
  expect_equal(res$richness, big$richness)
})

test_that("NULL mins/maxs fall back to band quantiles", {
  skip_without_gdal()
  skip_on_cran()

  set.seed(3)
  r <- terra::rast(nrows = 30, ncols = 30,
                   vals = runif(900, -2, 2))
  f <- tempfile(fileext = ".tif")
  terra::writeRaster(r, f, overwrite = TRUE)
  on.exit(unlink(f), add = TRUE)

  occ <- data.frame(cell = 1:30, species = "s1")
  qs <- terra::global(r, fun = function(v)
    stats::quantile(v, probs = c(0.005, 0.995), na.rm = TRUE))
  a <- pam_env_engine(f, occ, nbs = 4L)              # defaults path
  b <- pam_env_engine(f, occ, nbs = 4L,
                      mins = qs[[1]], maxs = qs[[2]])
  expect_equal(a$richness, b$richness)
  expect_equal(a$nnz, b$nnz)
})

test_that("has_gdal returns a logical flag", {
  expect_type(has_gdal(), "logical")
})

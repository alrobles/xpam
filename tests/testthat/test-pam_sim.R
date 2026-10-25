test_that("sim generator is deterministic", {
  skip_on_cran()
  r1 <- xpam:::xpam_sim_env(20, 15, 2, seed = 7)
  r2 <- xpam:::xpam_sim_env(20, 15, 2, seed = 7)
  expect_identical(terra::values(r1), terra::values(r2))
  s1 <- xpam:::xpam_sim_ranges(8, r1, seed = 3, mode = "range")
  s2 <- xpam:::xpam_sim_ranges(8, r1, seed = 3, mode = "range")
  expect_identical(s1$occ, s2$occ)
  s3 <- xpam:::xpam_sim_ranges(8, r1, seed = 4, mode = "range")
  expect_false(identical(s1$occ, s3$occ))
})

test_that("engine output equals transparent ground truth (niche mode)", {
  skip_without_gdal()
  skip_on_cran()

  nr <- 60; nc <- 45; D <- 2L; nbs <- c(5L, 6L)
  r <- xpam:::xpam_sim_env(nr, nc, D, seed = 1)
  f <- tempfile(fileext = ".tif")
  terra::writeRaster(r, f, overwrite = TRUE)
  on.exit(unlink(f), add = TRUE)

  sim <- xpam:::xpam_sim_ranges(15, r, seed = 2, mode = "niche")
  V <- terra::values(r, mat = TRUE)
  mins <- apply(V, 2, min) - 1e-9; maxs <- apply(V, 2, max) + 1e-9

  eng <- pam_env_engine(f, sim$occ, mins = mins, maxs = maxs, nbs = nbs,
                        tile_rows = 13, threads = 2)
  lv <- unique(sim$occ$species)
  truth <- xpam:::xpam_truth(V, sim$occ, lv, mins, maxs, nbs)

  expect_equal(eng$nnz, truth$nnz)
  expect_equal(eng$richness, truth$richness)
  expect_equal(eng$footprint, truth$footprint)
  expect_equal(eng$availability, truth$availability)

  # Analytic niche bound: footprint <= prod(per-axis hit bins)
  bound <- xpam:::xpam_niche_bound(sim$pars, mins, maxs, nbs)
  expect_true(all(eng$footprint <= bound[match(lv, sim$pars$species)]))
})

test_that("engine output equals transparent ground truth (range mode)", {
  skip_without_gdal()
  skip_on_cran()

  nr <- 50; nc <- 40; nbs <- c(4L, 4L)
  r <- xpam:::xpam_sim_env(nr, nc, 2, seed = 5)
  f1 <- tempfile(fileext = ".tif"); f2 <- tempfile(fileext = ".tif")
  terra::writeRaster(r[[1]], f1, overwrite = TRUE)
  terra::writeRaster(r[[2]], f2, overwrite = TRUE)
  on.exit(unlink(c(f1, f2)), add = TRUE)

  sim <- xpam:::xpam_sim_ranges(10, r, seed = 6, mode = "range")
  V <- terra::values(r, mat = TRUE)
  mins <- c(-4, -4); maxs <- c(4, 4)

  eng <- pam_env_engine(c(f1, f2), sim$occ, mins = mins, maxs = maxs,
                        nbs = nbs, tile_rows = 9, threads = 2)
  lv <- unique(sim$occ$species)
  truth <- xpam:::xpam_truth(V, sim$occ, lv, mins, maxs, nbs)

  expect_equal(eng$nnz, truth$nnz)
  expect_equal(eng$richness, truth$richness)
  expect_equal(eng$footprint, truth$footprint)
})

test_that("streamed occ_file equals in-memory occ (Phase C)", {
  skip_without_gdal()
  skip_on_cran()

  r <- xpam:::xpam_sim_env(80, 60, 2, seed = 11)
  f <- tempfile(fileext = ".tif")
  terra::writeRaster(r, f, overwrite = TRUE)
  on.exit(unlink(f), add = TRUE)

  sim <- xpam:::xpam_sim_ranges(25, r, seed = 12, mode = "range")
  V <- terra::values(r, mat = TRUE)
  mins <- apply(V, 2, min) - 1e-9; maxs <- apply(V, 2, max) + 1e-9
  nbs <- c(6L, 7L)

  a <- pam_env_engine(f, occ = sim$occ, mins = mins, maxs = maxs,
                      nbs = nbs, tile_rows = 11, threads = 2)

  fb <- tempfile(fileext = ".xpoc")
  on.exit(unlink(fb), add = TRUE)
  xpam:::xpam_write_occ_bin(sim$occ, fb)
  b <- pam_env_engine(f, occ_file = fb, mins = mins, maxs = maxs,
                      nbs = nbs, tile_rows = 11, threads = 2)

  for (k in c("nnz", "richness", "footprint", "availability", "n_strips"))
    expect_equal(b[[k]], a[[k]], label = paste("occ_file", k))
  expect_null(b$species)         # file carries codes, not names
})

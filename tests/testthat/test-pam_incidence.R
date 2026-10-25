# Sparse incidence export (issues #1/#2): the PAM itself as an R object.

test_that("geo incidence returns a units x species dgCMatrix", {
  skip_without_gdal()
  skip_on_cran()
  skip_if_not_installed("Matrix")

  nr <- 4; nc <- 5
  r <- terra::rast(nrows = nr, ncols = nc, vals = runif(nr * nc))
  f <- tempfile(fileext = ".tif")
  terra::writeRaster(r, f, overwrite = TRUE)
  on.exit(unlink(f), add = TRUE)

  # sp1: cells 1,2,3 | sp2: 2,3,4 | sp3: 5,6 (+ a duplicate pair)
  occ <- data.frame(
    cell = c(1, 2, 3, 2, 3, 4, 5, 6, 1),
    species = c("sp1", "sp1", "sp1", "sp2", "sp2", "sp2", "sp3", "sp3",
                "sp1"))

  res <- pam_incidence(f, occ, space = "geo")
  expect_s4_class(res$pam, "dgCMatrix")
  expect_equal(dim(res$pam), c(nr * nc, 3))
  expect_equal(colnames(res$pam), res$species)
  expect_equal(res$kernel, "bitmap")
  expect_equal(res$space, 1L)
  expect_equal(res$nnz, 8)              # duplicate (sp1, cell 1) counted once

  # incidence as a logical matrix
  M <- as.matrix(res$pam != 0)
  expect_equal(unname(M[, "sp1"]), 1:20 %in% c(1, 2, 3))
  expect_equal(unname(M[, "sp2"]), 1:20 %in% c(2, 3, 4))
  expect_equal(unname(M[, "sp3"]), 1:20 %in% c(5, 6))
  # marginals agree with the matrix
  expect_equal(unname(Matrix::colSums(res$pam)), res$footprint)
  expect_equal(sum(res$richness > 0), length(unique(occ$cell)))
})

test_that("format='ij' returns the raw 1-based pair table", {
  skip_without_gdal()
  skip_on_cran()

  nr <- 4; nc <- 5
  r <- terra::rast(nrows = nr, ncols = nc, vals = runif(nr * nc))
  f <- tempfile(fileext = ".tif")
  terra::writeRaster(r, f, overwrite = TRUE)
  on.exit(unlink(f), add = TRUE)

  occ <- data.frame(
    cell = c(1, 2, 3, 2, 3, 4, 5, 6, 1),
    species = c("sp1", "sp1", "sp1", "sp2", "sp2", "sp2", "sp3", "sp3",
                "sp1"))

  res <- pam_incidence(f, occ, space = "geo", format = "ij")
  expect_s3_class(res$pam, "data.frame")
  expect_equal(sort(colnames(res$pam)), c("i", "j"))
  sp_idx <- match(occ$species, res$species)
  got <- res$pam[order(res$pam$i, res$pam$j), ]
  want <- unique(data.frame(i = occ$cell, j = sp_idx))
  want <- want[order(want$i, want$j), ]
  expect_equal(unname(as.matrix(got)), unname(as.matrix(want)))
})

test_that("env-space incidence matches the binned pair table (sort path too)", {
  skip_without_gdal()
  skip_on_cran()
  skip_if_not_installed("Matrix")

  set.seed(7)
  r <- terra::rast(nrows = 10, ncols = 10, vals = runif(100))
  f <- tempfile(fileext = ".tif")
  terra::writeRaster(r, f, overwrite = TRUE)
  on.exit(unlink(f), add = TRUE)

  # 30 cells from 2 species; env space with D = 1, 10 bins over [0,1]
  occ <- data.frame(
    cell = sample(100, 30),
    species = rep(c("a", "b"), each = 15))
  vals <- terra::values(r)[, 1]
  bin_of <- function(cell)
    pmin(10L, pmax(1L, floor(vals[cell] * 10) + 1L))
  want <- unique(data.frame(i = bin_of(occ$cell),
                            j = match(occ$species, c("a", "b"))))
  want <- want[order(want$i, want$j), ]

  for (mb in c(1024, 0)) {                       # bitmap then sort fallback
    res <- pam_incidence(f, occ, space = "env",
                         mins = 0, maxs = 1, nbs = 10,
                         max_bitset_mb = mb)
    expect_equal(dim(res$pam), c(10, 2))
    M <- as.matrix(res$pam != 0)
    for (k in seq_len(nrow(want)))
      expect_true(M[want$i[k], want$j[k]])
    expect_equal(sum(M), nrow(want))
  }
})

test_that("geo over-budget refuses: co-occurrence retains no incidence", {
  skip_without_gdal()
  skip_on_cran()

  nr <- 4; nc <- 5
  r <- terra::rast(nrows = nr, ncols = nc, vals = runif(nr * nc))
  f <- tempfile(fileext = ".tif")
  terra::writeRaster(r, f, overwrite = TRUE)
  on.exit(unlink(f), add = TRUE)

  occ <- data.frame(cell = c(1, 2), species = c("a", "b"))
  expect_error(
    pam_incidence(f, occ, space = "geo", max_bitset_mb = 0),
    "co-occurrence|cooc|bitset|occ_dump")
})

test_that("nnz exceeding max_out_mb refuses the R return", {
  skip_without_gdal()
  skip_on_cran()

  nr <- 4; nc <- 5
  r <- terra::rast(nrows = nr, ncols = nc, vals = runif(nr * nc))
  f <- tempfile(fileext = ".tif")
  terra::writeRaster(r, f, overwrite = TRUE)
  on.exit(unlink(f), add = TRUE)

  occ <- data.frame(cell = 1:5, species = letters[1:5])
  expect_error(
    pam_incidence(f, occ, space = "geo", max_out_mb = 0),
    "max_out_mb|occ_dump|triplets")
})

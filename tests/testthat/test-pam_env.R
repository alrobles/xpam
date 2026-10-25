test_that("bitmap and sort kernels agree", {
  set.seed(1)
  ns <- 50; nb <- 7; P <- 2000
  sp <- sample.int(ns, P, TRUE) - 1L
  bn <- sample.int(nb, P, TRUE) - 1L
  a <- xpam:::c_env_pam_bitmap(sp, bn, ns, nb)
  b <- xpam:::c_env_pam(sp, bn, ns, nb)
  expect_equal(a$richness, b$richness)
  expect_equal(a$footprint, b$footprint)
  expect_equal(a$nnz, b$nnz)
})

test_that("bitmap kernel handles NAs and out-of-domain codes", {
  sp <- c(0L, 1L, NA_integer_, 2L)
  bn <- c(0L, NA_integer_, 1L, 1L)
  r <- xpam:::c_env_pam_bitmap(sp, bn, 3L, 2L)
  expect_equal(r$nnz, 2L)
  expect_equal(r$richness, c(1L, 1L))
  expect_equal(r$footprint, c(1L, 0L, 1L))
})

test_that("pam_env reproduces a manual bin count", {
  occ <- data.frame(cell = c(1, 1, 2, 3, 3, 4), species = c("a", "b", "a", "a", "c", "b"))
  env <- data.frame(cell = 1:4, x = c(0, 0, 1, 1), y = c(0, 1, 0, 1))
  r <- pam_env(occ, env, nb = 2, xlim = c(0, 1), ylim = c(0, 1))
  expect_equal(r$richness, c(2L, 1L, 2L, 1L))
  expect_equal(sort(r$footprint), c(1L, 2L, 3L))
  expect_equal(r$availability, c(1L, 1L, 1L, 1L))
})

test_that("shared-bitmap kernel loses no update under heavy word contention", {
  # The kernel ORs into ONE shared bitset from all threads. Concurrent writes to the
  # same 64-bit word are the failure mode (lost updates), so drive it with a grid whose
  # words are hit by many different (species, bin) pairs at once and cross-check against
  # the sort kernel, which is a single-threaded reference.
  set.seed(99)
  ns <- 4L; nbins <- 64L; P <- 5e5L
  sp <- sample.int(ns, P, TRUE) - 1L
  bn <- sample.int(nbins, P, TRUE) - 1L
  expect_equal(nbins, 64L)          # exactly one word per species row: worst case
  a <- xpam:::c_env_pam_bitmap(sp, bn, ns, nbins)
  b <- xpam:::c_env_pam(sp, bn, ns, nbins)
  expect_equal(a$nnz, b$nnz)
  expect_equal(a$richness, b$richness)
  expect_equal(a$footprint, b$footprint)
  expect_equal(a$n_out_of_range, 0L)
})

test_that("out-of-range codes are counted, not written, by the bitmap kernel", {
  r <- xpam:::c_env_pam_bitmap(c(0L, 5L, 1L), c(0L, 0L, 99L), 3L, 2L)
  expect_equal(r$n_out_of_range, 2L)
  expect_equal(r$nnz, 1L)
  expect_equal(r$richness, c(1L, 0L))
})

test_that("c_bin_assign_d generalises c_bin_assign to D dimensions", {
  set.seed(7)
  n <- 2000
  x <- runif(n, -3, 3); y <- rnorm(n)
  b2 <- xpam:::c_bin_assign(x, y, -3, 3, -4, 4, 7L)
  bd <- xpam:::c_bin_assign_d(cbind(x, y), c(-3, -4), c(3, 4), c(7L, 7L))
  expect_identical(b2, bd)                       # dim 0 most significant
  X3 <- cbind(x, y, rep(0:1, each = n / 2))
  b3 <- xpam:::c_bin_assign_d(X3, c(-3, -4, -1), c(3, 4, 2), c(7L, 7L, 2L))
  expect_true(all(b3 >= 0L & b3 < 7L * 7L * 2L, na.rm = TRUE))
  expect_true(all(is.na(xpam:::c_bin_assign_d(
    matrix(NA_real_, 2, 3), c(0, 0, 0), c(1, 1, 1), c(4L, 4L, 4L)))))
})

test_that("pam_env falls back to the sort kernel when the bitset exceeds the budget", {
  occ <- data.frame(cell = 1:6, species = rep(c("a", "b"), 3))
  env <- data.frame(cell = 1:6, x = seq(0, 1, length.out = 6), y = seq(0, 1, length.out = 6))
  r <- pam_env(occ, env, nb = 4, xlim = c(0, 1), ylim = c(0, 1), max_bitset_mb = 0)
  expect_equal(r$kernel, "sort")
  b <- pam_env(occ, env, nb = 4, xlim = c(0, 1), ylim = c(0, 1))
  expect_equal(b$kernel, "bitmap")
  expect_equal(r$richness, b$richness)
  expect_equal(r$footprint, b$footprint)
  expect_equal(r$nnz, b$nnz)
})

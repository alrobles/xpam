# R reference for pairwise Jaccard (fastJaccard::jaccard_r semantics):
# J = |A n B| / |A u B|; empty union -> expected pa*pb/(pa+pb-pa*pb);
# centered subtracts the expectation (empty union -> 0).
.jaccard_ref <- function(sets, U, centered = FALSE) {
  n <- length(sets)
  m <- matrix(NA_real_, n, n)
  for (i in seq_len(n)) {
    for (j in seq_len(i)) {
      a <- sets[[i]]; b <- sets[[j]]
      inter <- length(intersect(a, b))
      uni <- length(a) + length(b) - inter
      pa <- length(a) / U; pb <- length(b) / U
      e <- pa * pb / (pa + pb - pa * pb)
      v <- if (uni == 0) e else inter / uni
      if (centered) v <- if (uni == 0) 0 else v - e
      m[i, j] <- m[j, i] <- v
    }
  }
  m
}

test_that("geo-space Jaccard matches an R reference (bitset path)", {
  skip_without_gdal()
  skip_on_cran()

  nr <- 4; nc <- 5                                  # 20 cells
  r <- terra::rast(nrows = nr, ncols = nc, vals = runif(nr * nc))
  f <- tempfile(fileext = ".tif")
  terra::writeRaster(r, f, overwrite = TRUE)
  on.exit(unlink(f), add = TRUE)

  # sp1: cells 1,2,3 | sp2: 2,3,4 (J = 2/4) | sp3: 5,6 (disjoint -> 0)
  occ <- data.frame(
    cell = c(1, 2, 3, 2, 3, 4, 5, 6),
    species = c("sp1", "sp1", "sp1", "sp2", "sp2", "sp2", "sp3", "sp3"))

  res <- pam_jaccard(f, occ, space = "geo")
  expect_equal(res$kernel, "bitmap")
  expect_equal(res$space, 1L)
  expect_equal(dim(res$jaccard), c(3, 3))
  expect_equal(dimnames(res$jaccard)[[1]], res$species)

  sets <- split(occ$cell - 1L, factor(occ$species, levels = res$species))
  ref <- .jaccard_ref(sets, U = nr * nc)
  expect_equal(unname(res$jaccard), unname(ref), tolerance = 1e-12)
  expect_equal(res$jaccard[1, 2], 0.5, tolerance = 1e-12)
  expect_equal(res$jaccard[1, 3], 0)
  expect_equal(res$jaccard[3, 3], 1)
  # footprint = geographic range size in cells
  expect_equal(res$footprint, c(3L, 3L, 2L))
})

test_that("geo over-budget uses the co-occurrence accumulator", {
  skip_without_gdal()
  skip_on_cran()

  nr <- 4; nc <- 5
  r <- terra::rast(nrows = nr, ncols = nc, vals = runif(nr * nc))
  f <- tempfile(fileext = ".tif")
  terra::writeRaster(r, f, overwrite = TRUE)
  on.exit(unlink(f), add = TRUE)

  occ <- data.frame(
    cell = c(1, 2, 3, 2, 3, 4, 5, 6),
    species = c("sp1", "sp1", "sp1", "sp2", "sp2", "sp2", "sp3", "sp3"))

  res <- pam_jaccard(f, occ, space = "geo", max_bitset_mb = 0)
  expect_equal(res$kernel, "cooc")
  expect_equal(dim(res$cooc), c(3, 3))
  expect_equal(unname(diag(res$cooc)), c(3, 3, 2)) # diagonal = own cell count
  expect_equal(res$cooc[1, 2], 2)               # shared cells {2,3}
  expect_equal(res$cooc[1, 3], 0)
  expect_equal(res$jaccard[1, 2], 0.5, tolerance = 1e-12)
  expect_equal(res$jaccard[3, 3], 1)
})

test_that("env-space Jaccard agrees across bitmap and sort kernels", {
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

  jb <- pam_jaccard(c(f1, f2), occ, space = "env",
                    nbs = c(4, 5), mins = c(-4, -5), maxs = c(4, 5))
  expect_equal(jb$kernel, "bitmap")
  js <- pam_jaccard(c(f1, f2), occ, space = "env",
                    nbs = c(4, 5), mins = c(-4, -5), maxs = c(4, 5),
                    max_bitset_mb = 0)
  expect_equal(js$kernel, "sort")
  expect_equal(jb$jaccard, js$jaccard, tolerance = 1e-12)

  # R reference on the same 4x5 env grid
  bins <- xpam:::c_bin_assign_d(cbind(terra::values(r1)[, 1],
                                    terra::values(r2)[, 1]),
                              c(-4, -5), c(4, 5), c(4L, 5L))
  idx <- match(occ$cell, seq_len(nr * nc))
  bp <- bins[idx]
  ok <- bp >= 0
  sp <- factor(occ$species, levels = jb$species)
  sets <- lapply(levels(sp), function(s) sort(unique(bp[sp == s & ok])))
  ref <- .jaccard_ref(sets, U = 20)
  expect_equal(unname(jb$jaccard), unname(ref), tolerance = 1e-12)
})

test_that("centered mode subtracts the independence expectation", {
  skip_without_gdal()
  skip_on_cran()

  r <- terra::rast(nrows = 4, ncols = 5, vals = runif(20))
  f <- tempfile(fileext = ".tif")
  terra::writeRaster(r, f, overwrite = TRUE)
  on.exit(unlink(f), add = TRUE)

  occ <- data.frame(
    cell = c(1, 2, 3, 2, 3, 4, 5, 6),
    species = c("a", "a", "a", "b", "b", "b", "c", "c"))

  raw <- pam_jaccard(f, occ, space = "geo")
  cen <- pam_jaccard(f, occ, space = "geo", centered = TRUE)
  U <- 20
  pa <- pb <- 3 / U; pc <- 2 / U
  e_ab <- (pa * pb) / (pa + pb - pa * pb)
  expect_equal(cen$jaccard[1, 2], raw$jaccard[1, 2] - e_ab,
               tolerance = 1e-12)
  # non-identical correction per pair (pc differs)
  expect_false(isTRUE(all.equal(cen$jaccard[1, 2] - raw$jaccard[1, 2],
                                cen$jaccard[1, 3] - raw$jaccard[1, 3],
                                tolerance = 1e-12)))
})

test_that("duplicate pairs do not inflate co-occurrence counts", {
  skip_without_gdal()
  skip_on_cran()

  r <- terra::rast(nrows = 2, ncols = 3, vals = runif(6))
  f <- tempfile(fileext = ".tif")
  terra::writeRaster(r, f, overwrite = TRUE)
  on.exit(unlink(f), add = TRUE)

  occ <- data.frame(cell = c(1, 1, 1, 2, 2),
                    species = c("a", "a", "b", "a", "b"))  # dup (1,a)
  res <- pam_jaccard(f, occ, space = "geo", max_bitset_mb = 0)
  expect_equal(res$cooc[1, 2], 2)   # cells {1,2} shared once each
  expect_equal(unname(diag(res$cooc)), c(2, 2))
  expect_equal(res$jaccard[1, 2], 1, tolerance = 1e-12)
})

test_that("site-mode Jaccard over geo cells matches an R reference", {
  skip_without_gdal()
  skip_on_cran()

  nr <- 4; nc <- 5                                  # 20 units
  r <- terra::rast(nrows = nr, ncols = nc, vals = runif(nr * nc))
  f <- tempfile(fileext = ".tif")
  terra::writeRaster(r, f, overwrite = TRUE)
  on.exit(unlink(f), add = TRUE)

  # species per cell: c1{a}, c2{a,b}, c3{a,b}, c4{b}, c5{c}, c6{c}
  occ <- data.frame(
    cell = c(1, 2, 3, 2, 3, 4, 5, 6),
    species = c("a", "a", "a", "b", "b", "b", "c", "c"))

  res <- pam_jaccard(f, occ, space = "geo", by = "site")
  expect_equal(res$by, 1L)
  expect_equal(res$kernel, "bitmap")
  expect_equal(dim(res$jaccard), c(20, 20))         # units x units
  expect_equal(dimnames(res$jaccard)[[1]], as.character(1:20))

  sets <- lapply(seq_len(nr * nc), function(cell)
    sort(unique(occ$species[occ$cell == cell])))
  ref <- .jaccard_ref(sets, U = 3)                  # U = n_species
  expect_equal(unname(res$jaccard), unname(ref), tolerance = 1e-12)
  expect_equal(res$jaccard[2, 3], 1)                # identical {a,b}
  expect_equal(res$jaccard[1, 2], 0.5)
  expect_equal(res$jaccard[5, 6], 1)                # identical {c}
  expect_equal(res$jaccard[1, 5], 0)                # disjoint
  expect_true(is.nan(res$jaccard[7, 8]))            # both empty -> NaN

  # marginals keep public semantics under the transpose:
  # richness = species per unit, footprint = units per species
  expect_equal(res$richness[1:6], c(1L, 2L, 2L, 1L, 1L, 1L))
  expect_equal(res$footprint, c(3L, 3L, 2L))
})

test_that("site-mode Jaccard works over env bins", {
  skip_without_gdal()
  skip_on_cran()

  nr <- 20; nc <- 20
  set.seed(3)
  r1 <- terra::rast(nrows = nr, ncols = nc, vals = runif(nr * nc, 0, 1))
  r2 <- terra::rast(nrows = nr, ncols = nc, vals = runif(nr * nc, 0, 1))
  f1 <- tempfile(fileext = ".tif"); f2 <- tempfile(fileext = ".tif")
  terra::writeRaster(r1, f1, overwrite = TRUE)
  terra::writeRaster(r2, f2, overwrite = TRUE)
  on.exit(unlink(c(f1, f2)), add = TRUE)

  occ <- data.frame(
    cell = sample.int(nr * nc, 200, replace = TRUE),
    species = sample(paste0("sp", 1:8), 200, replace = TRUE))

  res <- pam_jaccard(c(f1, f2), occ, space = "env", by = "site",
                     nbs = c(3, 4), mins = c(0, 0), maxs = c(1, 1))
  expect_equal(res$by, 1L)
  expect_equal(dim(res$jaccard), c(12, 12))         # n_bins x n_bins

  bins <- xpam:::c_bin_assign_d(cbind(terra::values(r1)[, 1],
                                    terra::values(r2)[, 1]),
                              c(0, 0), c(1, 1), c(3L, 4L))
  idx <- match(occ$cell, seq_len(nr * nc))
  bp <- bins[idx]; ok <- bp >= 0
  sets <- lapply(0:11, function(b)
    sort(unique(occ$species[bp == b & ok])))
  ref <- .jaccard_ref(sets, U = 8)
  expect_equal(unname(res$jaccard), unname(ref), tolerance = 1e-12)
})

test_that("site mode fails cleanly when over budget", {
  skip_without_gdal()
  skip_on_cran()

  r <- terra::rast(nrows = 4, ncols = 5, vals = runif(20))
  f <- tempfile(fileext = ".tif")
  terra::writeRaster(r, f, overwrite = TRUE)
  on.exit(unlink(f), add = TRUE)

  occ <- data.frame(cell = c(1, 2, 2),
                    species = c("a", "a", "b"))

  # transposed bitset does not fit -> no fallback available
  expect_error(pam_jaccard(f, occ, space = "geo", by = "site",
                           max_bitset_mb = 0),
               "transposed bitset")
  # bitset fits but the 20x20 dense output (3.2 kB) exceeds max_out_mb
  expect_error(pam_jaccard(f, occ, space = "geo", by = "site",
                           max_out_mb = 1e-3),
               "dense Jaccard output")
})

# Collect the >= tau pairs of a dense symmetric matrix as (i, j, v),
# i < j, sorted — the oracle the edge list is compared against.
.edges_from_dense <- function(m, tau) {
  n <- nrow(m)
  out <- data.frame(i = integer(0), j = integer(0), jaccard = numeric(0))
  for (i in seq_len(n)) for (j in seq_len(i - 1))
    if (!is.na(m[i, j]) && m[i, j] >= tau)
      out <- rbind(out, data.frame(i = j, j = i, jaccard = m[i, j]))
  out[order(out$i, out$j), ]
}

test_that("threshold edges equal the dense matrix filtered at tau (geo)", {
  skip_without_gdal()
  skip_on_cran()

  nr <- 8; nc <- 10                                 # 80 cells
  set.seed(5)
  r <- terra::rast(nrows = nr, ncols = nc, vals = runif(nr * nc))
  f <- tempfile(fileext = ".tif")
  terra::writeRaster(r, f, overwrite = TRUE)
  on.exit(unlink(f), add = TRUE)

  # clustered species: a few share most of their cells
  occ <- data.frame(
    cell = c(1:20, 5:24, 10:29, 60:70),
    species = rep(paste0("sp", 1:4),
                  c(20, 20, 20, 11)))

  dense <- pam_jaccard(f, occ, space = "geo")
  for (tau in c(0.2, 0.5, 0.75, 1.0)) {
    res <- pam_jaccard(f, occ, space = "geo", threshold = tau)
    expect_null(res[["jaccard"]])
    expect_equal(res$threshold, tau)
    ref <- .edges_from_dense(unname(dense$jaccard), tau)
    expect_equal(res$edges$i, ref$i)
    expect_equal(res$edges$j, ref$j)
    expect_equal(res$edges$jaccard, ref$jaccard, tolerance = 1e-12)
  }
})

test_that("threshold edges work on sort and cooc accumulators", {
  skip_without_gdal()
  skip_on_cran()

  nr <- 8; nc <- 10
  set.seed(5)
  r <- terra::rast(nrows = nr, ncols = nc, vals = runif(nr * nc))
  f <- tempfile(fileext = ".tif")
  terra::writeRaster(r, f, overwrite = TRUE)
  on.exit(unlink(f), add = TRUE)

  occ <- data.frame(
    cell = c(1:20, 5:24, 10:29, 60:70),
    species = rep(paste0("sp", 1:4), c(20, 20, 20, 11)))
  tau <- 0.4

  dense <- pam_jaccard(f, occ, space = "geo")
  ref <- .edges_from_dense(unname(dense$jaccard), tau)

  # env space -> sort kernel when the bitset does not fit
  e_sort <- pam_jaccard(f, occ, space = "env", nbs = 4,
                        mins = 0, maxs = 1, threshold = tau,
                        max_bitset_mb = 0)
  expect_equal(e_sort$kernel, "sort")
  d_sort <- pam_jaccard(f, occ, space = "env", nbs = 4,
                      mins = 0, maxs = 1, max_bitset_mb = 0)
  ref_sort <- .edges_from_dense(unname(d_sort$jaccard), tau)
  expect_equal(e_sort$edges$i, ref_sort$i)
  expect_equal(e_sort$edges$jaccard, ref_sort$jaccard, tolerance = 1e-12)

  # geo over-budget -> cooc kernel
  e_cooc <- pam_jaccard(f, occ, space = "geo", threshold = tau,
                        max_bitset_mb = 0)
  expect_equal(e_cooc$kernel, "cooc")
  expect_equal(e_cooc$edges$i, ref$i)
  expect_equal(e_cooc$edges$j, ref$j)
  expect_equal(e_cooc$edges$jaccard, ref$jaccard, tolerance = 1e-12)
})

test_that("threshold edges work in site mode and when empty", {
  skip_without_gdal()
  skip_on_cran()

  nr <- 4; nc <- 5
  r <- terra::rast(nrows = nr, ncols = nc, vals = runif(nr * nc))
  f <- tempfile(fileext = ".tif")
  terra::writeRaster(r, f, overwrite = TRUE)
  on.exit(unlink(f), add = TRUE)

  occ <- data.frame(
    cell = c(1, 2, 3, 2, 3, 4, 5, 6),
    species = c("a", "a", "a", "b", "b", "b", "c", "c"))

  res <- pam_jaccard(f, occ, space = "geo", by = "site", threshold = 0.9)
  expect_s3_class(res$edges, "data.frame")
  # cells 2,3 hold {a,b}; cells 5,6 hold {c} -> exactly two J=1 edges
  expect_equal(res$edges$i, c(2, 5))
  expect_equal(res$edges$j, c(3, 6))
  expect_equal(res$edges$jaccard, c(1, 1))
  # tau = 1 keeps the same edges (identical species sets)
  res1 <- pam_jaccard(f, occ, space = "geo", by = "site", threshold = 1)
  expect_equal(res1$edges$i, c(2, 5))

  # every cell holds a different species -> no edges at all
  occ0 <- data.frame(cell = 1:3, species = c("a", "b", "c"))
  res0 <- pam_jaccard(f, occ0, space = "geo", by = "site", threshold = 1)
  expect_null(res0[["jaccard"]])
  expect_equal(nrow(res0$edges), 0)
})

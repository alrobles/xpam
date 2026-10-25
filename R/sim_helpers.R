# sim_helpers.R -- deterministic synthetic worlds for tests and benchmarks.
#
# xpam_sim_env:    smooth deterministic env fields -> SpatRaster (D layers)
# xpam_sim_ranges: species (cell, species) pairs; mode "niche" = boxes in env
#                  space (analytic footprint bound), mode "range" = geographic
#                  ellipses (IUCN-like).
# xpam_truth:      transparent ground-truth recomputation in plain R.
#
# All functions are internal (not exported) and fully deterministic given
# their seed.

#' @keywords internal
#' @noRd
xpam_sim_env <- function(nr, nc, D = 2L, seed = 1L) {
  x <- rep(seq(-1, 1, length.out = nc), times = nr)
  y <- rep(seq(-1, 1, length.out = nr), each = nc)
  set.seed(seed)
  vals <- matrix(NA_real_, nr * nc, D)
  for (d in seq_len(D)) {
    z <- switch((d - 1L) %% 4L + 1L,
      x,                                   # d 1: lon gradient
      y,                                   # d 2: lat gradient
      x * y,                               # d 3: interaction
      sin(3 * x) * cos(3 * y))             # d 4: patchy wave
    if (d > 4L) z <- z + runif(nr * nc, -0.5, 0.5)   # extra dims get noise
    vals[, d] <- as.numeric(scale(z))
  }
  r <- terra::rast(nrows = nr, ncols = nc, nlyrs = D,
                   extent = c(0, nc, 0, nr))
  for (d in seq_len(D)) terra::values(r[[d]]) <- vals[, d]
  names(r) <- paste0("env", seq_len(D))
  r
}

#' @keywords internal
#' @noRd
xpam_sim_ranges <- function(n_species, env_rast, seed = 1L,
                            mode = c("niche", "range"),
                            tol = c(0.25, 0.8), area = c(0.005, 0.05)) {
  mode <- match.arg(mode)
  set.seed(seed)
  nr <- terra::nrow(env_rast); nc <- terra::ncol(env_rast)
  if (mode == "niche") {
    V <- terra::values(env_rast, mat = TRUE)
    mu  <- matrix(runif(n_species * 2, -0.9, 0.9), n_species, 2)
    tl  <- matrix(runif(n_species * 2, tol[1], tol[2]), n_species, 2)
    cells <- vector("list", n_species)
    for (s in seq_len(n_species)) {
      inside <- abs(V[, 1] - mu[s, 1]) <= tl[s, 1] &
                abs(V[, 2] - mu[s, 2]) <= tl[s, 2]
      cells[[s]] <- which(inside)
    }
    pars <- data.frame(species = paste0("sp", seq_len(n_species)),
                       mu1 = mu[, 1], mu2 = mu[, 2],
                       tol1 = tl[, 1], tol2 = tl[, 2])
  } else {
    cx <- runif(n_species, 0.1 * nc, 0.9 * nc)
    cy <- runif(n_species, 0.1 * nr, 0.9 * nr)
    ar <- runif(n_species, area[1], area[2]) * nr * nc
    rx <- sqrt(ar / pi); ry <- rx * runif(n_species, 0.5, 2)
    xx <- rep(seq_len(nc), times = nr); yy <- rep(seq_len(nr), each = nc)
    cells <- vector("list", n_species)
    for (s in seq_len(n_species)) {
      inside <- ((xx - cx[s]) / rx[s])^2 + ((yy - cy[s]) / ry[s])^2 <= 1
      cells[[s]] <- which(inside)
    }
    pars <- data.frame(species = paste0("sp", seq_len(n_species)),
                       cx = cx, cy = cy, rx = rx, ry = ry)
  }
  occ <- data.frame(
    cell = unlist(cells),
    species = rep(pars$species, lengths(cells)))
  list(occ = occ, pars = pars)
}

#' Analytic footprint bound for mode = "niche": the niche box
#' [mu_d - tol_d, mu_d + tol_d] intersects at most B_d bins on axis d, so
#' footprint <= prod_d B_d (equality when every intersected bin column is
#' populated).  Exact recomputation still happens via xpam_truth().
#' @keywords internal
#' @noRd
xpam_niche_bound <- function(pars, mins, maxs, nbs) {
  hit <- function(lo, hi, mn, mx, nb) {
    i <- pmin(pmax(floor((lo - mn) / (mx - mn) * nb), 0), nb - 1)
    j <- pmin(pmax(floor((hi - mn) / (mx - mn) * nb), 0), nb - 1)
    j - i + 1
  }
  b1 <- hit(pars$mu1 - pars$tol1, pars$mu1 + pars$tol1, mins[1], maxs[1], nbs[1])
  b2 <- hit(pars$mu2 - pars$tol2, pars$mu2 + pars$tol2, mins[2], maxs[2], nbs[2])
  b1 * b2
}

#' Transparent ground truth: bin every cell in plain R (no C++), deduplicate
#' (species, bin) pairs, return marginals.  Independent of every kernel.
#' @keywords internal
#' @noRd
xpam_truth <- function(env_vals, occ, species_lv, mins, maxs, nbs) {
  D <- length(nbs)
  stride <- rev(cumprod(rev(nbs))) / nbs         # stride_d = prod_{k>d} nb_k
  bin_of <- function(row) {
    id <- 0
    for (d in seq_len(D)) {
      v <- row[d]
      if (!is.finite(v)) return(-1L)
      k <- floor((v - mins[d]) / (maxs[d] - mins[d]) * nbs[d])
      id <- id + min(max(k, 0), nbs[d] - 1) * stride[d]
    }
    id
  }
  nb_total <- prod(nbs)
  cellbin <- apply(env_vals, 1L, bin_of)
  rich <- integer(nb_total); foot <- integer(length(species_lv))
  sp_code <- match(occ$species, species_lv)
  seen <- new.env(hash = TRUE, parent = emptyenv())
  nnz <- 0L
  for (p in seq_len(nrow(occ))) {
    b <- cellbin[occ$cell[p]]
    if (b < 0L) next
    key <- paste(sp_code[p], b, sep = ":")
    if (is.null(seen[[key]])) {
      seen[[key]] <- TRUE; nnz <- nnz + 1L
      rich[b + 1L] <- rich[b + 1L] + 1L
      foot[sp_code[p]] <- foot[sp_code[p]] + 1L
    }
  }
  avail <- integer(nb_total)
  for (b in cellbin) if (b >= 0L) avail[b + 1L] <- avail[b + 1L] + 1L
  list(richness = rich, footprint = foot, availability = as.numeric(avail),
       nnz = nnz)
}

#' Write an occurrence table to the streamed "XPOC" v1 binary format the
#' engine reads with FileOccSource: header + float64 cell-sorted cells +
#' int32 species codes.  Cells are terra 1-based cellnumbers here (the R
#' conductor converts to 0-based), so we subtract 1 while writing.
#' @keywords internal
#' @noRd
xpam_write_occ_bin <- function(occ, path, species_lv = NULL) {
  if (is.null(species_lv)) species_lv <- unique(occ$species)
  sp_code <- as.integer(match(occ$species, species_lv)) - 1L
  cells <- as.double(occ$cell) - 1.0                    # terra -> 0-based
  keep <- !is.na(cells) & !is.na(sp_code)
  cells <- cells[keep]; sp_code <- sp_code[keep]
  o <- order(cells)
  con <- file(path, "wb"); on.exit(close(con))
  writeChar("XPOC", con, eos = NULL)
  writeBin(1L, con, size = 4L)                          # version
  writeBin(c(length(cells), length(species_lv)), con, size = 8L)
  writeBin(cells[o], con, size = 8L)                    # float64, sorted
  writeBin(sp_code[o], con, size = 4L)                  # int32
  invisible(species_lv)
}

#!/usr/bin/env Rscript
# bench_engine_scaling.R -- one scaling configuration per invocation.
#
# Usage:
#   Rscript bench_engine_scaling.R <side> <n_species> <threads> <mode> [cover]
#
#   side       raster side in pixels (grid is side x side)
#   n_species  species in the simulated pool
#   threads    OpenMP threads for the engine
#   mode       "prep" | "engine" | "vilela"
#   cover      range-size span as fraction of side (default 0.125); shrink it
#              at large sides or the pair table grows as cells x species
#
# prep   : build deterministic env rasters + occurrence table once, cache to
#          bench_cache/<side>_<nsp>/ (untimed; dominated by R-side work).
# engine : read the cached inputs, run pam_env_engine, print one csv row.
# vilela : dense-paradigm baseline — materialize the sites x species PAM and
#          map to env bins (the lets.presab approach).  Reported while the
#          dense double matrix fits; beyond that it prints NA — which IS the
#          point of the comparison.
#
# Each run is a fresh process so VmHWM is a clean per-config peak RSS
# (Linux /proc/self/status) measuring the algorithm, not data prep.

suppressMessages({library(terra); library(xpam)})
args <- commandArgs(trailingOnly = TRUE)
side <- as.integer(args[1]); nsp <- as.integer(args[2])
nth <- as.integer(args[3]); mode <- args[4]
cover <- if (length(args) >= 5) as.numeric(args[5]) else 0.03
gate_gb <- if (length(args) >= 6) as.numeric(args[6]) else 3
nr <- nc <- side; D <- 2L; nbs <- c(30L, 30L)

dir.create(cache <- file.path("bench_cache",
                              sprintf("%d_%d_%.3f", side, nsp, cover)),
           recursive = TRUE, showWarnings = FALSE)
fenv <- file.path(cache, "env.tif"); focc <- file.path(cache, "occ.rds")

rss_mb <- function() {
  v <- grep("VmHWM", readLines("/proc/self/status"), value = TRUE)
  as.numeric(gsub("[^0-9]", "", v)) / 1024
}
catrow <- function(mode, secs, rss, nnz = NA_real_) {
  cat(sprintf("%s,%d,%d,%d,%d,%.3f,%.1f,%s\n", mode, side,
              as.double(nr) * nc, nsp, nth, secs, rss,
              ifelse(is.na(nnz), "NA", format(nnz))))
}

if (mode == "prep") {
  r <- xpam:::xpam_sim_env(nr, nc, D, seed = 1)
  writeRaster(r, fenv, overwrite = TRUE)
  sim <- xpam:::xpam_sim_ranges(nsp, r, seed = 2, mode = "range",
                               area = c(0.3, 1.7) * cover)
  V <- values(r, mat = TRUE)
  mins <- apply(V, 2, min) - 1e-9; maxs <- apply(V, 2, max) + 1e-9
  saveRDS(list(occ = sim$occ, pars = sim$pars), focc)
  saveRDS(list(mins = mins, maxs = maxs, pars = sim$pars),
          file.path(cache, "meta.rds"))
  xpam:::xpam_write_occ_bin(sim$occ, file.path(cache, "occ.xpoc"))
  cat("prep done", cache, "\n")
  quit("no")
}

meta <- readRDS(file.path(cache, "meta.rds"))
mins <- meta$mins; maxs <- meta$maxs; cfg <- meta

if (mode == "engine") {
  occ <- readRDS(focc)$occ
  t0 <- proc.time()[3]
  res <- pam_env_engine(fenv, occ = occ, mins = mins, maxs = maxs, nbs = nbs,
                        tile_rows = 512, threads = nth)
  catrow("engine", proc.time()[3] - t0, rss_mb(), res$nnz)
} else if (mode == "engine_file") {
  # pair table stays on disk entirely -- engine streams it
  t0 <- proc.time()[3]
  res <- pam_env_engine(fenv, occ_file = file.path(cache, "occ.xpoc"),
                        mins = mins, maxs = maxs, nbs = nbs,
                        tile_rows = 512, threads = nth)
  catrow("engine_file", proc.time()[3] - t0, rss_mb(), res$nnz)
} else {
  dense_mb <- as.double(nr) * nc * nsp * 8 / 1024^2   # real double matrix
  if (dense_mb > gate_gb * 1024) {
    catrow("vilela", NA_real_, NA_real_)
  } else {
    t0 <- proc.time()[3]
    occ <- readRDS(focc)$occ
    lay <- vector("list", nsp)
    for (s in seq_len(nsp)) {
      cells <- occ$cell[occ$species == cfg$pars$species[s]]
      lr <- terra::rast(nrows = nr, ncols = nc, vals = 0L)
      lr[cells] <- 1L
      lay[[s]] <- lr
    }
    pam_dense <- rast(lay)
    V <- values(rast(fenv), mat = TRUE)
    bins <- xpam:::c_bin_assign_d(V, mins, maxs, nbs)
    ok <- bins >= 0L
    M <- values(pam_dense, mat = TRUE)[ok, , drop = FALSE]
    B <- matrix(0L, sum(ok), prod(nbs))
    B[cbind(seq_len(nrow(B)), bins[ok] + 1L)] <- 1L
    epam <- crossprod(M, B) > 0
    catrow("vilela", proc.time()[3] - t0, rss_mb(), sum(epam))
  }
}

#!/usr/bin/env Rscript
# sweep_run.R -- xpam real-data sweep, one configuration per invocation.
#
# Usage:
#   Rscript sweep_run.R <res> <mode> <threads> [year]
#
#   res      "t10"  (1/10 deg, coarse50 bio01+bio12) |
#            "p120" (1/120 deg, era5land 1 km bio01+bio12)
#   mode     "poly" (IUCN shapefile -> PolygonOccSource, dumps XPOC) |
#            "file" (stream the dumped XPOC via occ_file) |
#            "mem"  (read the XPOC into an R data.frame -> occ)
#   threads  OpenMP threads for the engine
#   year     bioclim year (default 2020)
#
# Each run is a fresh process so VmHWM is a clean per-config peak RSS.
# Prints one CSV row: mode,res,n_cells,n_species,threads,seconds,rss_mb,
# nnz,n_skipped.  Run "poly" before "file"/"mem": it writes the pair
# table the other modes stream.
#
# Paths are cluster-local (/beegfs); sweep.sh wraps the config grid.

suppressMessages({library(terra); library(xpam)})
args <- commandArgs(trailingOnly = TRUE)
res_id <- args[1]; mode <- args[2]; nth <- as.integer(args[3])
year <- if (length(args) >= 4) as.integer(args[4]) else 2020L

ROOT  <- "/beegfs/a474r867/xpam_sweep"
WORK  <- file.path(ROOT, "work")
DIV   <- "/beegfs/a474r867/diversity_work"
STAGE <- "/beegfs/a474r867/stage1_global_hybrid/bioclim_separate_1km_global_1986-2025"

env_files <- switch(res_id,
  t10  = file.path(DIV, "coarse50",
                   sprintf("bio%02d_%d.tif", c(1, 12), year)),
  p120 = file.path(STAGE, as.character(year),
                   sprintf("bio%02d_%d.tif", c(1, 12), year)),
  conus = file.path(WORK, sprintf("conus_bio%02d_%d.tif", c(1, 12), year)),
  stop("unknown res: ", res_id))
stopifnot(all(file.exists(env_files)))
gpkg <- file.path(WORK, "mammals.gpkg")
xpoc <- Sys.getenv("XPOC",
                   file.path(WORK, sprintf("occ_%s.xpoc", res_id)))
meta <- readRDS(file.path(WORK, sprintf("meta_%s.rds", res_id)))
nbs  <- c(30L, 30L)                       # bio01 x bio12, paper grid

rss_mb <- function() {
  v <- grep("VmHWM", readLines("/proc/self/status"), value = TRUE)
  as.numeric(gsub("[^0-9]", "", v)) / 1024
}
catrow <- function(mode, secs, rss, res) {
  nsp <- if (length(res$species)) length(res$species)
         else length(res$footprint)   # one footprint entry per species
  cat(sprintf("%s,%s,%s,%d,%d,%.3f,%.1f,%s,%s\n", mode, res_id,
              format(meta$ncells, scientific = FALSE), nsp, nth, secs, rss,
              ifelse(is.null(res$nnz), "NA", format(res$nnz)),
              ifelse(is.null(res$n_skipped), "NA",
                     format(res$n_skipped))))
}

if (mode == "poly") {
  ## The pair stream is identical across thread counts -- only the
  ## first poly run pays the dump write (at p120 it is ~100 GB).
  t0 <- proc.time()[3]
  res <- pam_env_engine(env_files, occ_polygons = gpkg,
                        poly_field = "sci_name", all_touched = TRUE,
                        occ_dump = if (nth == 1L) xpoc else NULL,
                        mins = meta$mins, maxs = meta$maxs, nbs = nbs,
                        tile_rows = 512, threads = nth)
  catrow("poly", proc.time()[3] - t0, rss_mb(), res)
} else if (mode == "file") {
  t0 <- proc.time()[3]
  res <- pam_env_engine(env_files, occ_file = xpoc,
                        mins = meta$mins, maxs = meta$maxs, nbs = nbs,
                        tile_rows = 512, threads = nth)
  catrow("file", proc.time()[3] - t0, rss_mb(), res)
} else if (mode == "mem") {
  con <- file(xpoc, "rb")
  stopifnot(readChar(con, 4, useBytes = TRUE) == "XPOC")
  readBin(con, "integer", 1, size = 4)
  ## np/nsp are little-endian int64: assemble from raw bytes so values
  ## beyond 2^31 are not silently truncated to the low word.
  rd64 <- function() sum(as.integer(readBin(con, "raw", 8L)) * 256^(0:7))
  np  <- rd64()
  nsp <- rd64()
  ## Gate BEFORE materialising: a resident pair table costs ~40 B/row
  ## in R (float64 + int32 vectors plus data.frame row.names and GC
  ## copies -- 9.2e8 pairs OOM-killed a 24 G job at conus; 12 B/pair
  ## underestimates), and data.frames are capped at 2^31-1 rows anyway.
  ## Beyond ~24 GB the run cannot fit the allocation -> report
  ## infeasible (the NA row IS the result) instead of OOM-dying mid-read.
  est_gb <- np * 40 / 1e9
  if (np >= 2^31 - 1 || est_gb > 24) {
    message(sprintf("mem mode infeasible: np=%.3g pairs, ~%.1f GB",
                    np, est_gb))
    cat(sprintf("mem,%s,%s,%s,%d,NA,NA,NA,NA\n", res_id,
                format(meta$ncells, scientific = FALSE),
                format(nsp, scientific = FALSE), nth))
    quit(save = "no", status = 0)
  }
  occ <- data.frame(cell = readBin(con, "double", np, size = 8) + 1,
                    species = readBin(con, "integer", np, size = 4))
  close(con)
  t0 <- proc.time()[3]
  res <- pam_env_engine(env_files, occ = occ, mins = meta$mins,
                        maxs = meta$maxs, nbs = nbs,
                        tile_rows = 512, threads = nth)
  res$species <- seq_len(nsp)
  catrow("mem", proc.time()[3] - t0, rss_mb(), res)
} else {
  stop("unknown mode: ", mode)
}

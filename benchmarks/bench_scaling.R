# Scaling benchmark: bitmap kernel vs sort-based kernel vs naive R (data.table).
# Usage: Rscript benchmarks/bench_scaling.R [data_dir] [max_reps]
suppressMessages({library(data.table); library(xpam)})
args <- commandArgs(trailingOnly = TRUE)
dd <- if (length(args) >= 1) args[1] else "/home/alrobles/GitHub/diversity_assesment/data/iucn"
maxrep <- if (length(args) >= 2) as.integer(args[2]) else 16L

cells <- fread(file.path(dd, "iucn_cells.csv")); occ <- fread(file.path(dd, "iucn_occurrence.csv"))
ser <- fread(file.path(dd, "iucn_clim_series.csv"))[year == 2000][is.finite(bio01) & is.finite(bio12)]
cel <- merge(cells, ser[, .(cell_id, bio01, bio12)], by = "cell_id")
z1 <- (cel$bio01 - mean(cel$bio01)) / sd(cel$bio01); z2 <- (cel$bio12 - mean(cel$bio12)) / sd(cel$bio12)
xlim <- range(z1); ylim <- range(z2)
NB <- 30; n_bins <- NB * NB
n_species <- uniqueN(occ$species)
bins <- xpam:::c_bin_assign(z1, z2, xlim[1], xlim[2], ylim[1], ylim[2], NB)
celR <- data.table(cell_id = cel$cell_id, bid = bins)[!is.na(bid)]

tm <- function(expr) { g <- gc(reset = TRUE); t <- system.time(force(expr))[["elapsed"]]; m <- sum(gc()[, 6]); c(t = t, mem = m) }

cat(sprintf("n_species %d  n_cells %d  n_bins %d\n", n_species, nrow(cel), n_bins))
cat(sprintf("dense species x cell matrix if materialised: %.1f GB  (the reason fine grids are infeasible)\n",
            n_species * nrow(cel) * 8 / 1e9))
cat(sprintf("bitmap species x bin state: %.2f MB  (n_species*n_bins/8)\n\n", n_species * n_bins / 8 / 1e6))

cat(sprintf("%10s %12s %12s %12s %14s %12s\n", "reps", "pairs", "R_dt_s", "bitmap_s", "sort_s", "bitmap_MB"))
for (reps in c(1L, 4L, 8L, maxrep)) {
  o <- occ[rep(seq_len(.N), reps)]
  sp <- match(o$species, unique(o$species)) - 1L
  bp <- bins[match(o$cell_id, cel$cell_id)]
  a <- tm(unique(merge(o, celR, by = "cell_id")[, .(bid, species)]))
  b <- tm(xpam:::c_env_pam_bitmap(sp, bp, n_species, n_bins))
  s <- tm(xpam:::c_env_pam(sp, bp, n_species, n_bins))
  cat(sprintf("%10d %12d %12.2f %12.2f %12.2f %14.1f\n", reps, nrow(o), a["t"], b["t"], s["t"], b["mem"]))
}
cat("\n(R_dt_s = naive data.table merge+unique; bitmap/sort = C++ kernels; bitmap_MB = peak R memory in the bitmap run)\n")

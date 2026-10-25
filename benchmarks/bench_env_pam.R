# Benchmark: sparse C++ environmental PAM vs the naive R (merge + unique) path.
# Usage: Rscript benchmarks/bench_env_pam.R <data_dir>
# <data_dir> holds iucn_cells.csv, iucn_occurrence.csv, iucn_clim_series.csv.

suppressMessages({library(data.table); library(xpam)})
args <- commandArgs(trailingOnly = TRUE)
dd <- if (length(args) >= 1) args[1] else "/home/alrobles/GitHub/diversity_assesment/data/iucn"

cells <- fread(file.path(dd, "iucn_cells.csv"))
occ   <- fread(file.path(dd, "iucn_occurrence.csv"))
ser   <- fread(file.path(dd, "iucn_clim_series.csv"))[year == 2000]
ser   <- ser[is.finite(bio01) & is.finite(bio12)]
cel   <- merge(cells, ser[, .(cell_id, bio01, bio12)], by = "cell_id")

z1 <- (cel$bio01 - mean(cel$bio01)) / sd(cel$bio01)
z2 <- (cel$bio12 - mean(cel$bio12)) / sd(cel$bio12)
xlim <- range(z1); ylim <- range(z2); NB <- 30

cat(sprintf("cells %d  pairs %d  species %d\n", nrow(cel), nrow(occ), uniqueN(occ$species)))

## ---- naive R : merge + unique + count ---------------------------------------
t0 <- Sys.time()
bins <- xpam:::c_bin_assign(z1, z2, xlim[1], xlim[2], ylim[1], ylim[2], NB)   # same binning as C++
celR <- data.table(cell_id = cel$cell_id, bid = bins)[!is.na(bid)]
sb   <- unique(merge(occ, celR, by = "cell_id")[, .(bid, species)])
rich_R <- sb[, .N, by = bid]
fp_R   <- sb[, .N, by = species]
t_R <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

## ---- C++ sparse path ---------------------------------------------------------
t0 <- Sys.time()
res <- pam_env(occ = data.frame(cell = occ$cell_id, species = occ$species),
               env = data.frame(cell = cel$cell_id, x = z1, y = z2),
               nb = NB, xlim = xlim, ylim = ylim)
t_C <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

## ---- compare ----------------------------------------------------------------
rich_C <- data.table(bid = res$j, N = 1)[, .N, by = bid][order(bid)]
setorder(rich_R, bid)
ok_rich <- nrow(rich_R) == nrow(rich_C) && all(rich_R$bid == rich_C$bid) && all(rich_R$N == rich_C$N)
cat(sprintf("\nrichness identical (bin-level): %s   bins %d vs %d\n",
            ok_rich, nrow(rich_R), nrow(rich_C)))
cat(sprintf("footprint sizes identical: %s\n",
            identical(sort(fp_R$N), sort(res$footprint))))
cat(sprintf("nnz: R %d  C++ %d\n", nrow(sb), res$nnz))
cat(sprintf("\ntime  naive R: %.2f s   C++: %.2f s   speedup %.1fx\n", t_R, t_C, t_R / t_C))

## ---- scaling probe: replicate pairs to emulate a finer grid -----------------
occ2 <- occ[rep(seq_len(.N), 4)]
cat(sprintf("\nscaling probe: %d pairs\n", nrow(occ2)))
t0 <- Sys.time()
invisible(unique(merge(occ2, celR, by = "cell_id")[, .(bid, species)]))
t_R2 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
t0 <- Sys.time()
invisible(pam_env(occ = data.frame(cell = occ2$cell_id, species = occ2$species),
                  env = data.frame(cell = cel$cell_id, x = z1, y = z2),
                  nb = NB, xlim = xlim, ylim = ylim))
t_C2 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("  naive R: %.2f s   C++: %.2f s   speedup %.1fx\n", t_R2, t_C2, t_R2 / t_C2))

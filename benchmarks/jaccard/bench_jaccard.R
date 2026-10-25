# bench_jaccard.R -- pairwise Jaccard kernel benchmark.
#
# Compares on the SAME synthetic species x cell incidence:
#   1. fastJaccard::jaccard_fast_matrix   dense doubles, RcppParallel
#   2. xpam::pam_jaccard (dense)          bitset rows, OpenMP
#   3. xpam::pam_jaccard (threshold)      exact prefix-filtered join
#
# The fixture mimics ecological incidence: most species occupy small,
# spatially coherent ranges over a large grid, a few are widespread --
# the regime where a threshold join should dominate the dense O(n^2)
# pass.
#
# Usage:
#   R_LIBS=/tmp/xlib Rscript benchmarks/jaccard/bench_jaccard.R
# (xpam and fastJaccard installed in R_LIBS; results append to
# benchmarks/jaccard/bench_results.csv)

suppressPackageStartupMessages({
  library(xpam)
  library(terra)
})

has_fj <- requireNamespace("fastJaccard", quietly = TRUE)

set.seed(7)
n_sp <- 2000L          # entities
nr <- 100L; nc <- 200L # U = 20,000 cells
U <- nr * nc

# --- synthetic incidence -------------------------------------------------
# species ranges = elliptical blobs around a random centre + rare
# widespread species; sizes log-normal ~ realistic range-size spectrum.
mk_range <- function() {
  if (runif(1) < 0.02) {                       # ~2% widespread
    return(sample.int(U, size = sample(2000:8000, 1)))
  }
  cx <- runif(1, 1, nc); cy <- runif(1, 1, nr)
  a <- runif(1, 2, 18); b <- runif(1, 2, 18)   # semi-axes
  xx <- rep(1:nc, times = nr); yy <- rep(1:nr, each = nc)
  d2 <- ((xx - cx) / a)^2 + ((yy - cy) / b)^2
  cells <- which(d2 <= 1)
  if (length(cells) == 0) cells <- as.integer(cy * nc + cx) + 1L
  cells
}
sets <- lapply(seq_len(n_sp), function(i) mk_range())

# template raster (geometry only) + occ pairs for the engine path
tpl <- tempfile(fileext = ".tif")
writeRaster(rast(nrows = nr, ncols = nc, vals = seq_len(U)), tpl,
            overwrite = TRUE)
occ <- data.frame(
  cell    = unlist(sets),
  species = rep(paste0("sp", seq_len(n_sp)), lengths(sets)))
cat(sprintf("fixture: %d species x %d cells, %d pairs\n",
            n_sp, U, nrow(occ)))

# dense binary matrix for fastJaccard (rows = species)
mat <- matrix(0.0, nrow = n_sp, ncol = U)
for (s in seq_len(n_sp)) mat[s, sets[[s]]] <- 1
cat(sprintf("dense input matrix: %.1f MB\n",
            object.size(mat) / 1048576))

rows <- list()
stamp <- function(algo, sec, extra) {
  rows[[length(rows) + 1]] <<- data.frame(
    algo = algo, n_entities = n_sp, n_units = U,
    seconds = round(sec, 3), extra = extra,
    stringsAsFactors = FALSE)
  cat(sprintf("%-38s %8.2f s   %s\n", algo, sec, extra))
}

# --- 1. fastJaccard dense (baseline to beat) ------------------------------
if (has_fj) {
  t <- system.time(fjm <- fastJaccard::jaccard_fast_matrix(mat))[["elapsed"]]
  stamp("fastJaccard::jaccard_fast_matrix", t,
        sprintf("dist matrix %dx%d", nrow(fjm), ncol(fjm)))
  rm(fjm); gc()
} else {
  cat("fastJaccard not installed - skipping baseline\n")
}
rm(mat); gc()

# --- 2. xpam dense (bitset rows) ------------------------------------------
t <- system.time(res_d <- pam_jaccard(
  tpl, occ, space = "geo", threads = max(1L, parallel::detectCores() %/% 2L)
))[["elapsed"]]
stamp("xpam dense (bitset, e2e)", t,
      sprintf("kernel=%s matrix=%dx%d", res_d$kernel,
              nrow(res_d$jaccard), ncol(res_d$jaccard)))

# --- 3. xpam thresholded exact join ---------------------------------------
for (tau in c(0.05, 0.2, 0.5)) {
  t <- system.time(res_e <- pam_jaccard(
    tpl, occ, space = "geo", threshold = tau,
    threads = max(1L, parallel::detectCores() %/% 2L)))[["elapsed"]]
  stamp(sprintf("xpam edges tau=%.2f (e2e)", tau), t,
        sprintf("kernel=%s edges=%d", res_e$kernel, nrow(res_e$edges)))
}

# --- cross-check: edges at tau equal dense filtered -----------------------
tau_chk <- 0.2
res_c <- pam_jaccard(tpl, occ, space = "geo", threshold = tau_chk,
                     threads = max(1L, parallel::detectCores() %/% 2L))
m <- unname(res_d$jaccard)
keep <- upper.tri(m) & (m >= tau_chk)          # canonical i < j, like edges
ii <- row(m)[keep]; jj <- col(m)[keep]; vv <- m[keep]
o <- order(ii, jj)
stopifnot(identical(res_c$edges$i, ii[o]),
          identical(res_c$edges$j, jj[o]),
          isTRUE(all.equal(res_c$edges$jaccard, vv[o],
                           tolerance = 1e-12)))
cat("edge-vs-dense cross-check at tau=0.2: OK\n")

out <- do.call(rbind, rows)
csv <- "benchmarks/jaccard/bench_results.csv"
if (file.exists(csv)) {
  old <- read.csv(csv)
  out <- rbind(old, out[!(out$algo %in% old$algo & out$n_entities == n_sp), ])
}
write.csv(out, csv, row.names = FALSE)
cat("wrote", csv, "\n")

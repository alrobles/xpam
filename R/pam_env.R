#' Build a sparse environmental presence-absence matrix
#'
#' Constructs the species-by-climate-bin incidence matrix (the "environmental
#' PAM") from a species occurrence table and a climate table, without ever
#' materialising a dense species-by-site or species-by-bin matrix. Climate space
#' is discretised onto an \code{nb x nb} grid; each cell is assigned to a bin and
#' each species is recorded as present in the bins its cells fall in.
#'
#' The accumulation is done in C++ with a single shared bitset
#' (\code{n_species * n_bins / 8} bytes, filled with lock-free atomic ORs), so
#' peak memory is independent of both the number of cells and the number of
#' threads, and no per-thread state has to be merged. This is what allows the
#' construction to run at resolutions (e.g. 1/120 degree) where the dense path is
#' infeasible.
#'
#' The kernel is chosen from a memory budget, and the choice is reported. If the
#' bitset for the requested binning would exceed \code{max_bitset_mb} MiB, the
#' sort-based kernel is used instead: it has a much smaller constant when
#' \code{n_species * n_bins} is huge, but its memory is \code{O(P)} (proportional
#' to the number of (species, cell) pairs), which is what the bitset avoids.
#' Raise the budget, or lower \code{nb}, if the bitset is affordable.
#'
#' @param occ A data frame with columns \code{cell} and \code{species}. \code{cell}
#'   must be matchable to \code{env$cell}; \code{species} identifies the taxon.
#' @param env A data frame with columns \code{cell}, \code{x}, \code{y}: the two
#'   climate coordinates of each cell.
#' @param nb Integer. Number of bins per climate axis. Default \code{30}.
#' @param xlim,ylim Numeric length-2 ranges for the binning. If \code{NULL},
#'   the 0.5\% and 99.5\% quantiles of the data are used.
#' @param max_bitset_mb Numeric. Memory budget (MiB) for the bitset kernel; if
#'   the bitset would exceed it, the sort-based kernel is used instead. Default
#'   \code{1024}.
#'
#' @return A list with \code{i}, \code{j} (0-based species and bin indices of the
#'   nonzero entries), \code{nnz}, \code{n_species}, \code{n_bins},
#'   \code{richness} (distinct species per bin), \code{footprint} (bins per
#'   species), \code{availability} (cells per bin, the climate-availability
#'   measure), \code{kernel} (the kernel actually used), and \code{state_mb}
#'   (the bitset size the binning would need).
#'
#' @export
#'
#' @examples
#' \donttest{
#' occ <- data.frame(cell = c(1, 1, 2, 3), species = c("a", "b", "a", "c"))
#' env <- data.frame(cell = 1:3, x = c(0.1, 0.4, 0.8), y = c(0.2, 0.5, 0.9))
#' pam_env(occ, env, nb = 4)
#' }
pam_env <- function(occ, env, nb = 30, xlim = NULL, ylim = NULL, max_bitset_mb = 1024) {
  checkmate::assert_data_frame(occ)
  checkmate::assert_data_frame(env)
  checkmate::assert_number(max_bitset_mb, lower = 0)
  stopifnot(all(c("cell", "species") %in% names(occ)))
  stopifnot(all(c("cell", "x", "y") %in% names(env)))
  nb <- as.integer(nb)
  stopifnot(nb >= 1L)

  if (is.null(xlim)) xlim <- as.numeric(stats::quantile(env$x, c(0.005, 0.995), na.rm = TRUE))
  if (is.null(ylim)) ylim <- as.numeric(stats::quantile(env$y, c(0.005, 0.995), na.rm = TRUE))

  bins <- c_bin_assign(env$x, env$y, xlim[1], xlim[2], ylim[1], ylim[2], nb)

  idx <- match(occ$cell, env$cell)
  bin_of_pair <- bins[idx]
  species_lv <- unique(occ$species)
  spp <- match(occ$species, species_lv) - 1L
  n_species <- length(species_lv)
  n_bins <- nb * nb

  # The bitset kernel costs n_species * n_bins / 8 bytes REGARDLESS of the thread
  # count; the previous per-thread-bitmap version multiplied that by the number of
  # threads, which is why the old guard compared against 2e9 bits (250 MiB per
  # thread) and silently fell back to the O(P) sort kernel above it.
  state_mb <- as.numeric(n_species) * n_bins / 8 / 1024^2
  use_bitmap <- state_mb <= max_bitset_mb

  if (use_bitmap) {
    res <- c_env_pam_bitmap(spp, bin_of_pair, n_species, n_bins)
    res$i <- NULL; res$j <- NULL
    if (res$n_out_of_range > 0L) {
      warning(sprintf("%d (species, bin) pair(s) fell outside the declared ranges and were skipped",
                      res$n_out_of_range))
    }
    res$n_out_of_range <- NULL
    kernel <- "bitmap"
  } else {
    res <- c_env_pam(spp, bin_of_pair, n_species, n_bins)
    res$i <- NULL; res$j <- NULL
    kernel <- "sort"
  }

  res$availability <- c_availability(bins, n_bins)
  res$species <- species_lv
  res$nb <- nb
  res$kernel <- kernel
  res$state_mb <- state_mb
  res
}

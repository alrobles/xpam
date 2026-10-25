// pam_kernel.cpp -- sparse, streaming environmental-PAM kernels for the pam package.
//
// The package's dense path (rasterize -> sites x species tibble) cannot scale to fine
// grids: a 1/120-deg global grid is ~9e8 cells and a species-by-site dense matrix is
// impossible. These kernels keep everything sparse: they accumulate only realised
// (species, climate-bin) incidences, so memory is proportional to the output, not to
// (cells x species).
//
// This file is a thin Rcpp layer: the heavy state lives in the pure-C++ headers
// (pam_bitset.h, pam_binning.h, pam_accumulate.h, pam_threads.h), which the
// engine reuses without Rcpp. All codes are 0-based integers (NA_INTEGER skipped).
//
// Thread counts go through pam::safe_threads(), which honours OMP_THREAD_LIMIT
// (omp_get_max_threads alone does not, and CRAN check machines set it to 2).

#include <Rcpp.h>
#include <vector>
#include <cstdint>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "pam_bitset.h"
#include "pam_binning.h"
#include "pam_accumulate.h"
#include "pam_threads.h"

using namespace Rcpp;

// Assign each point to a bin on an nb x nb grid spanning [xmin,xmax] x [ymin,ymax].
// Returns 0-based bin id = bx * nb + by, or NA for non-finite input.
// [[Rcpp::export]]
IntegerVector c_bin_assign(NumericVector x, NumericVector y,
                           double xmin, double xmax, double ymin, double ymax, int nb) {
  const R_xlen_t n = x.size();
  IntegerVector out(n);
  pam::GridBinning gb({xmin, ymin}, {xmax, ymax}, {nb, nb});
  const int nt = pam::safe_threads();
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nt)
#endif
  for (R_xlen_t k = 0; k < n; ++k) {
    const double xy[2] = {x[k], y[k]};
    const int64_t id = gb.assign(xy);
    out[k] = id < 0 ? NA_INTEGER : (int)id;
  }
  return out;
}

// D-dimensional version: X is n x D (column-major), mins/maxs/nbs length D.
// Dim 0 is the most significant digit, so for D = 2 this matches c_bin_assign.
// [[Rcpp::export]]
IntegerVector c_bin_assign_d(NumericMatrix X, NumericVector mins,
                             NumericVector maxs, IntegerVector nbs) {
  const R_xlen_t n = X.nrow();
  const int D = X.ncol();
  if (mins.size() != D || maxs.size() != D || nbs.size() != D)
    stop("mins, maxs and nbs must all have length ncol(X)");
  pam::GridBinning gb(as<std::vector<double>>(mins),
                      as<std::vector<double>>(maxs),
                      as<std::vector<int>>(nbs));
  IntegerVector out(n);
  const double* Xp = X.begin();
  const int nt = pam::safe_threads();
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nt)
#endif
  for (R_xlen_t k = 0; k < n; ++k) {
    const int64_t id = gb.assign_strided(Xp + k, n);
    out[k] = id < 0 ? NA_INTEGER : (int)id;
  }
  return out;
}

// Availability measure: number of cells per climate bin.
// [[Rcpp::export]]
IntegerVector c_availability(IntegerVector bin_of_cell, int n_bins) {
  IntegerVector a(n_bins);
  const R_xlen_t n = bin_of_cell.size();
  for (R_xlen_t k = 0; k < n; ++k) {
    const int b = bin_of_cell[k];
    if (b != NA_INTEGER && b >= 0 && b < n_bins) a[b]++;
  }
  return a;
}

// Sparse binary species x bin incidence (the environmental PAM), deduplicated,
// returned in coordinate form plus marginals. Keys are packed into uint64 and
// sorted; memory is O(number of input pairs), never O(n_species * n_bins).
// [[Rcpp::export]]
List c_env_pam(IntegerVector species, IntegerVector bin, int n_species, int n_bins) {
  const R_xlen_t P = species.size();
  if (bin.size() != P) stop("species and bin must have the same length");
  std::vector<uint64_t> keys;
  keys.reserve((size_t)P);
  for (R_xlen_t k = 0; k < P; ++k) {
    const int s = species[k], b = bin[k];
    if (s == NA_INTEGER || b == NA_INTEGER) continue;
    if (s < 0 || s >= n_species) stop("species code out of range");
    if (b < 0 || b >= n_bins) stop("bin code out of range");
    keys.push_back(((uint64_t)(uint32_t)s << 32) | (uint32_t)b);
  }
  std::sort(keys.begin(), keys.end());
  keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
  const R_xlen_t nnz = (R_xlen_t)keys.size();
  IntegerVector i(nnz), j(nnz);
  IntegerVector footprint(n_species);   // bins per species
  IntegerVector richness(n_bins);       // species per bin
  for (R_xlen_t k = 0; k < nnz; ++k) {
    const int s = (int)(keys[k] >> 32);
    const int b = (int)(keys[k] & 0xffffffffu);
    i[k] = s; j[k] = b; footprint[s]++; richness[b]++;
  }
  return List::create(_["i"] = i, _["j"] = j,
                      _["n_species"] = n_species, _["n_bins"] = n_bins,
                      _["nnz"] = nnz, _["footprint"] = footprint, _["richness"] = richness);
}

// Richness field: number of DISTINCT species per bin, given a long (species, cell)
// table and a cell->bin map. Streams in place, no dense intermediate.
// [[Rcpp::export]]
IntegerVector c_richness(IntegerVector species, IntegerVector bin_of_pair, int n_bins) {
  const R_xlen_t P = species.size();
  std::vector<uint64_t> keys;
  keys.reserve((size_t)P);
  for (R_xlen_t k = 0; k < P; ++k) {
    const int s = species[k], b = bin_of_pair[k];
    if (s == NA_INTEGER || b == NA_INTEGER) continue;
    keys.push_back(((uint64_t)(uint32_t)b << 32) | (uint32_t)s);
  }
  std::sort(keys.begin(), keys.end());
  keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
  IntegerVector rich(n_bins);
  for (size_t k = 0; k < keys.size(); ++k) rich[(int)(keys[k] >> 32)]++;
  return rich;
}

// Bitmap kernel: ONE shared species x bin bitset (pam::PamBitset), filled with
// lock-free atomic ORs via pam::accumulate_pairs, then marginals via population
// counts. Memory is n_species*n_bins/8 bytes TOTAL -- independent of the thread
// count and of the number of input pairs.
// [[Rcpp::export]]
List c_env_pam_bitmap(IntegerVector species, IntegerVector bin, int n_species, int n_bins) {
  const R_xlen_t P = species.size();
  if (bin.size() != P) stop("species and bin must have the same length");
  pam::PamBitset bs((size_t)n_species, (size_t)n_bins);
  const int64_t oor = pam::accumulate_pairs(bs, INTEGER(species), INTEGER(bin),
                                            (int64_t)P, NA_INTEGER,
                                            pam::safe_threads());
  IntegerVector rich(n_bins), foot(n_species);
  const size_t nnz = bs.marginals(INTEGER(rich), INTEGER(foot));
  return List::create(_["n_species"] = n_species, _["n_bins"] = n_bins,
                      _["nnz"] = (R_xlen_t)nnz, _["richness"] = rich, _["footprint"] = foot,
                      _["n_out_of_range"] = (R_xlen_t)oor);
}

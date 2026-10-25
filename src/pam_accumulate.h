// pam_accumulate.h -- the streaming scatter shared by the R kernel and the
// future engine: (species, bin) pair arrays -> PamBitset, with out-of-range
// pairs counted per thread instead of throwing (the R API must not be called
// from inside a parallel region). NA_INTEGER in either stream is skipped.
//
// Pure C++17 (no Rcpp); the only R-ism is the NA_INTEGER sentinel value
// itself (-2147483648), passed in explicitly so callers from non-R contexts
// can use their own sentinel.
#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>
#include "pam_bitset.h"
#include "pam_threads.h"

namespace pam {

// Accumulate P pairs of 0-based codes into bs. Returns the total number of
// pairs that were skipped for being out of range (NAs are skipped silently,
// matching the historical kernel behaviour).
inline std::int64_t accumulate_pairs(PamBitset& bs, const int* sp, const int* bn,
                                     std::int64_t P, int na_code,
                                     int nthreads = safe_threads()) {
  const std::int64_t S = (std::int64_t)bs.n_species();
  const std::int64_t B = (std::int64_t)bs.n_bins();
  std::int64_t oor = 0;
#ifdef _OPENMP
  if (nthreads > 1) {
    std::vector<std::int64_t> oort((std::size_t)nthreads, 0);
#pragma omp parallel num_threads(nthreads)
    {
      const int t = omp_get_thread_num();
#pragma omp for schedule(static)
      for (std::int64_t k = 0; k < P; ++k) {
        const int s = sp[k], b = bn[k];
        if (s == na_code || b == na_code) continue;
        if (s < 0 || s >= S || b < 0 || b >= B) { oort[(std::size_t)t]++; continue; }
        bs.add_atomic((std::size_t)s, (std::size_t)b);
      }
    }
    for (int t = 0; t < nthreads; ++t) oor += oort[(std::size_t)t];
    return oor;
  }
#endif
  (void)nthreads;
  for (std::int64_t k = 0; k < P; ++k) {
    const int s = sp[k], b = bn[k];
    if (s == na_code || b == na_code) continue;
    if (s < 0 || s >= S || b < 0 || b >= B) { oor++; continue; }
    bs.add((std::size_t)s, (std::size_t)b);
  }
  return oor;
}

} // namespace pam

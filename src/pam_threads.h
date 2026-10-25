// pam_threads.h -- thread-count helper honouring OMP_THREAD_LIMIT.
//
// omp_get_max_threads() respects OMP_NUM_THREADS but NOT OMP_THREAD_LIMIT;
// CRAN check machines set OMP_THREAD_LIMIT=2, and exceeding it during checks
// is a flagged failure. Same pattern as xbioclim_safe_threads(). Pure C++.
#pragma once

#include <cstdlib>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace pam {

// requested < 1 means "use the runtime maximum" (omp_get_max_threads()).
// The result is always capped by OMP_THREAD_LIMIT when that var is set.
inline int safe_threads(int requested = 0) noexcept {
#ifdef _OPENMP
  int n = requested > 0 ? requested : omp_get_max_threads();
  const char* lim = std::getenv("OMP_THREAD_LIMIT");
  if (lim && *lim) {
    const long v = std::strtol(lim, nullptr, 10);
    if (v > 0 && v < n) n = (int)v;
  }
  return n > 0 ? n : 1;
#else
  (void)requested;
  return 1;
#endif
}

} // namespace pam

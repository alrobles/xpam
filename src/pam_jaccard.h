// pam_jaccard.h -- pairwise Jaccard similarity inside the engine.
//
// Semantics ported from fastJaccard::jaccard_fast (alrobles/fastJaccard):
//   J(A,B) = |A n B| / |A u B|
//   union == 0  ->  E(pa,pb) = pa*pb / (pa + pb - pa*pb)  (expected value)
//   centered    ->  J - E(pa,pb)
// with pa = |A| / U, U = number of incidence units (env bins or cells).
//
// Implementation is new: packed bitset rows (popcount), sorted unit lists
// (merge join) or a co-occurrence count matrix -- never a dense double
// matrix, and OpenMP via pam::safe_threads() instead of RcppParallel.
//
// All kernels fill a full symmetric n x n matrix in `out` (row-major),
// except the *_edges variants, which emit only pairs meeting a threshold.
#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <vector>
#include "pam_bitset.h"
#include "pam_threads.h"

namespace pam {

// Probabilistic expectation of Jaccard under independence (the fallback
// value for empty unions and the centring term).  NaN when both sets are
// empty (pa = pb = 0).
inline double jaccard_expected(double pa, double pb) noexcept {
  const double den = pa + pb - pa * pb;
  return den > 0.0 ? (pa * pb) / den
                   : std::numeric_limits<double>::quiet_NaN();
}

// One pairwise value: raw J, expectation on union == 0, centred on demand.
inline double jaccard_value(double inter, double uni, double pa, double pb,
                            bool centered) noexcept {
  const double e = jaccard_expected(pa, pb);
  if (uni <= 0.0) return centered ? 0.0 : e;
  const double j = inter / uni;
  return centered ? j - e : j;
}

// Bitset rows: J[a,b] = popcount(A & B) / popcount(A | B).
inline void jaccard_rows(const PamBitset& bs, bool centered, int threads,
                         std::vector<double>& out) {
  const std::size_t n = bs.n_species();
  const std::size_t w = bs.words_per_row();
  const double U = (double)bs.n_bins();
  out.assign(n * n, 0.0);
  std::vector<std::uint64_t> sz(n, 0);
  for (std::size_t s = 0; s < n; ++s) {
    const std::uint64_t* r = bs.row(s);
    for (std::size_t k = 0; k < w; ++k) sz[s] += (std::uint64_t)__builtin_popcountll(r[k]);
  }
  const int nt = pam::safe_threads(threads);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4) num_threads(nt)
#endif
  for (std::int64_t i = 0; i < (std::int64_t)n; ++i) {
    const std::uint64_t* A = bs.row((std::size_t)i);
    const double pa = (double)sz[(std::size_t)i] / U;
    for (std::size_t j = 0; j <= (std::size_t)i; ++j) {
      const std::uint64_t* B = bs.row(j);
      double inter = 0.0, uni = 0.0;
      for (std::size_t k = 0; k < w; ++k) {
        inter += (double)__builtin_popcountll(A[k] & B[k]);
        uni += (double)__builtin_popcountll(A[k] | B[k]);
      }
      const double v = jaccard_value(inter, uni, pa, (double)sz[j] / U, centered);
      out[(std::size_t)i * n + j] = v;
      out[j * n + (std::size_t)i] = v;
    }
  }
}

// Sorted unit-id lists (sort-fallback path): merge join per pair.
inline void jaccard_lists(const std::vector<std::vector<std::int32_t>>& sets,
                          std::size_t n_units, bool centered, int threads,
                          std::vector<double>& out) {
  const std::size_t n = sets.size();
  const double U = (double)n_units;
  out.assign(n * n, 0.0);
  const int nt = pam::safe_threads(threads);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4) num_threads(nt)
#endif
  for (std::int64_t i = 0; i < (std::int64_t)n; ++i) {
    const std::vector<std::int32_t>& A = sets[(std::size_t)i];
    const double pa = (double)A.size() / U;
    for (std::size_t j = 0; j <= (std::size_t)i; ++j) {
      const std::vector<std::int32_t>& B = sets[j];
      std::size_t a = 0, b = 0, inter = 0;
      while (a < A.size() && b < B.size()) {
        if (A[a] == B[b]) { ++inter; ++a; ++b; }
        else if (A[a] < B[b]) ++a;
        else ++b;
      }
      const double uni = (double)(A.size() + B.size() - inter);
      const double v = jaccard_value((double)inter, uni, pa,
                                     (double)B.size() / U, centered);
      out[(std::size_t)i * n + j] = v;
      out[j * n + (std::size_t)i] = v;
    }
  }
}

// ── Thresholded edge output ─────────────────────────────────────────────
//
// Exact prefix-filtered similarity join (AllPairs-style, Bayardo et al.
// SIGMOD'07): for J >= tau the output is the edge list (i, j, v), not a
// dense matrix, and the work is output-sensitive instead of O(n^2).
//
// Derivation.  J(A,B) = o/(a+b-o) >= tau  <=>  o >= beta*(a+b) with
// beta = tau/(1+tau).  Entities are processed in increasing size order;
// entity y of size b is indexed only on the (b - t_y + 1) globally
// least-frequent units of its set, where t_y = ceil(2*beta*b) is the
// smallest overlap y could need with any partner >= its own size.  Its
// unindexed suffix then holds t_y - 1 units.  For a true match (x, y),
// b <= a:  o >= t = ceil(beta*(a+b)) >= t_y, so at least one shared unit
// lies in y's indexed prefix -- probing the inverted index with the full
// set of x cannot miss it.  Candidates are pruned by the size filter
// (b < tau*a cannot match) and the count bound (cnt + t_y - 1 < t), then
// verified by exact merge join -- every emitted edge is exact.
//
// For centred output the raw-J bound is still exact: v_centred >= tau
// implies v_raw >= tau (expectation >= 0), so generating candidates at
// raw tau and emitting only centred values >= tau loses no edge.
//
// sets: per-entity unit ids sorted ascending (used for the exact verify);
// freq: global unit counts (richness per unit) used to pick the prefix;
// empty sets never match tau > 0 and are skipped.
inline void jaccard_edges(
    const std::vector<std::vector<std::int32_t>>& sets,
    const std::vector<std::int64_t>& freq,
    std::size_t n_units, double tau, bool centered, int threads,
    std::vector<std::int32_t>& ei, std::vector<std::int32_t>& ej,
    std::vector<double>& ev) {
  const std::size_t n = sets.size();
  ei.clear(); ej.clear(); ev.clear();
  if (tau <= 0.0 || tau > 1.0 || n < 2) return;
  const double beta = tau / (1.0 + tau);
  // All required-overlap bounds below are integers in exact arithmetic;
  // subtract a tiny epsilon before ceil() so a mathematically integral
  // bound is not rounded up by last-ulp error (e.g. beta*276 = 46.0000..4
  // -> ceil 47 would wrongly drop pairs sitting exactly at J = tau).
  const double eps = 1e-6;

  // Entities in increasing set-size order; drop the empty ones.
  std::vector<std::size_t> order;
  order.reserve(n);
  for (std::size_t i = 0; i < n; ++i)
    if (!sets[i].empty()) order.push_back(i);
  std::sort(order.begin(), order.end(), [&sets](std::size_t p, std::size_t q) {
    return sets[p].size() < sets[q].size();
  });
  if (order.size() < 2) return;

  // Inverted index over indexed prefixes only (sparse: only units that
  // actually occur become keys).
  std::unordered_map<std::int32_t, std::vector<std::int32_t>> inv;

  std::vector<std::int32_t> cnt(n, 0), seen(n, 0), touched;
  touched.reserve(256);
  std::vector<std::pair<std::int32_t, std::int32_t>> cand;
  cand.reserve(1u << 16);

  std::vector<std::int32_t> byfreq;      // entity's units, rare-first

  for (std::size_t pos = 0; pos < order.size(); ++pos) {
    const std::size_t x = order[pos];
    const std::vector<std::int32_t>& xs = sets[x];
    const std::int64_t a = (std::int64_t)xs.size();
    const std::int64_t min_b =
      (std::int64_t)std::ceil(tau * (double)a - eps);

    // Probe: count, for every earlier entity y, |x n prefix(y)|.
    touched.clear();
    for (std::size_t k = 0; k < xs.size(); ++k) {
      auto it = inv.find(xs[k]);
      if (it == inv.end()) continue;
      for (std::int32_t y : it->second) {
        if (!seen[(std::size_t)y]) {
          seen[(std::size_t)y] = 1;
          touched.push_back(y);
          cnt[(std::size_t)y] = 1;
        } else {
          ++cnt[(std::size_t)y];
        }
      }
    }
    for (std::int32_t y : touched) {
      seen[(std::size_t)y] = 0;                       // reset for next x
      const std::int64_t b = (std::int64_t)sets[(std::size_t)y].size();
      if (b < min_b) { cnt[(std::size_t)y] = 0; continue; }   // size filter
      const std::int64_t t =
        (std::int64_t)std::ceil(beta * (double)(a + b) - eps);
      const std::int64_t t_y =
        (std::int64_t)std::ceil(2.0 * beta * (double)b - eps);
      if ((int64_t)cnt[(std::size_t)y] + t_y - 1 < t) {        // count bound
        cnt[(std::size_t)y] = 0;
        continue;
      }
      cand.emplace_back((std::int32_t)x, y);
      cnt[(std::size_t)y] = 0;
    }

    // Index x's least-frequent prefix: p(x) = a - t_x + 1 elements.
    const std::int64_t t_x =
      (std::int64_t)std::ceil(2.0 * beta * (double)a - eps);
    const std::size_t p = (std::size_t)(a - t_x + 1);
    if (p > 0) {
      byfreq.assign(xs.begin(), xs.end());
      const std::size_t pk = std::min(p, byfreq.size());
      std::partial_sort(byfreq.begin(), byfreq.begin() + (std::ptrdiff_t)pk,
                        byfreq.end(), [&freq](std::int32_t u, std::int32_t v) {
        const std::int64_t fu = u < (std::int32_t)freq.size() ? freq[(std::size_t)u] : 0;
        const std::int64_t fv = v < (std::int32_t)freq.size() ? freq[(std::size_t)v] : 0;
        return fu < fv || (fu == fv && u < v);
      });
      for (std::size_t k = 0; k < pk; ++k) inv[byfreq[k]].push_back((std::int32_t)x);
    }
  }

  // Verify candidates exactly and emit edges meeting the threshold
  // (parallel; per-thread buffers merged in candidate order after sort).
  const int nt = pam::safe_threads(threads);
  const double U = (double)n_units;
  std::vector<std::vector<std::int32_t>> ti((std::size_t)nt),
    tj((std::size_t)nt);
  std::vector<std::vector<double>> tv((std::size_t)nt);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 256) num_threads(nt)
#endif
  for (std::int64_t c = 0; c < (std::int64_t)cand.size(); ++c) {
    const std::size_t x = (std::size_t)cand[(std::size_t)c].first;
    const std::size_t y = (std::size_t)cand[(std::size_t)c].second;
    const std::vector<std::int32_t>& A = sets[x];
    const std::vector<std::int32_t>& B = sets[y];
    std::size_t ai = 0, bi = 0, inter = 0;
    while (ai < A.size() && bi < B.size()) {
      if (A[ai] == B[bi]) { ++inter; ++ai; ++bi; }
      else if (A[ai] < B[bi]) ++ai;
      else ++bi;
    }
    const double uni = (double)(A.size() + B.size() - inter);
    const double pa = (double)A.size() / U, pb = (double)B.size() / U;
    const double v = jaccard_value((double)inter, uni, pa, pb, centered);
    if (v >= tau) {
      const std::size_t th =
#ifdef _OPENMP
        (std::size_t)omp_get_thread_num();
#else
        0;
#endif
      ti[th].push_back((std::int32_t)std::min(x, y));
      tj[th].push_back((std::int32_t)std::max(x, y));
      tv[th].push_back(v);
    }
  }
  for (int th = 0; th < nt; ++th) {
    ei.insert(ei.end(), ti[(std::size_t)th].begin(), ti[(std::size_t)th].end());
    ej.insert(ej.end(), tj[(std::size_t)th].begin(), tj[(std::size_t)th].end());
    ev.insert(ev.end(), tv[(std::size_t)th].begin(), tv[(std::size_t)th].end());
  }
  // Deterministic edge order (i, j).
  std::vector<std::size_t> ord(ei.size());
  for (std::size_t k = 0; k < ord.size(); ++k) ord[k] = k;
  std::sort(ord.begin(), ord.end(), [&ei, &ej](std::size_t p, std::size_t q) {
    return ei[p] < ei[q] || (ei[p] == ei[q] && ej[p] < ej[q]);
  });
  std::vector<std::int32_t> si(ei.size()), sj(ej.size());
  std::vector<double> sv(ev.size());
  for (std::size_t k = 0; k < ord.size(); ++k) {
    si[k] = ei[ord[k]]; sj[k] = ej[ord[k]]; sv[k] = ev[ord[k]];
  }
  ei.swap(si); ej.swap(sj); ev.swap(sv);
}

// Thresholded edges straight from the co-occurrence accumulator: still an
// O(n^2) scan (that path is the memory fallback anyway) but the output
// stays sparse -- no dense matrix is materialised.
inline void jaccard_edges_from_cooc(
    const std::vector<std::uint64_t>& C,
    const std::vector<std::int64_t>& tot,
    std::size_t n, std::size_t n_units,
    double tau, bool centered, int threads,
    std::vector<std::int32_t>& ei, std::vector<std::int32_t>& ej,
    std::vector<double>& ev) {
  ei.clear(); ej.clear(); ev.clear();
  const double U = (double)n_units;
  const int nt = pam::safe_threads(threads);
  std::vector<std::vector<std::int32_t>> ti((std::size_t)nt), tj((std::size_t)nt);
  std::vector<std::vector<double>> tv((std::size_t)nt);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 32) num_threads(nt)
#endif
  for (std::int64_t i = 0; i < (std::int64_t)n; ++i) {
    const double pa = (double)tot[(std::size_t)i] / U;
    for (std::size_t j = 0; j < (std::size_t)i; ++j) {
      const double inter = (double)C[(std::size_t)i * n + j];
      const double uni = (double)tot[(std::size_t)i] + (double)tot[j] - inter;
      const double v = jaccard_value(inter, uni, pa, (double)tot[j] / U,
                                     centered);
      if (v >= tau) {
        const std::size_t th =
#ifdef _OPENMP
          (std::size_t)omp_get_thread_num();
#else
          0;
#endif
        ti[th].push_back((std::int32_t)j);
        tj[th].push_back((std::int32_t)i);
        tv[th].push_back(v);
      }
    }
  }
  for (int th = 0; th < nt; ++th) {
    ei.insert(ei.end(), ti[(std::size_t)th].begin(), ti[(std::size_t)th].end());
    ej.insert(ej.end(), tj[(std::size_t)th].begin(), tj[(std::size_t)th].end());
    ev.insert(ev.end(), tv[(std::size_t)th].begin(), tv[(std::size_t)th].end());
  }
}

// Co-occurrence counts: C[a*n+b] = units shared; tot[s] = units per species
// (geo path above the bitset budget; the C diagonal holds tot[a]).
inline void jaccard_from_cooc(const std::vector<std::uint64_t>& C,
                              const std::vector<std::int64_t>& tot,
                              std::size_t n, std::size_t n_units,
                              bool centered, int threads,
                              std::vector<double>& out) {
  const double U = (double)n_units;
  out.assign(n * n, 0.0);
  const int nt = pam::safe_threads(threads);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4) num_threads(nt)
#endif
  for (std::int64_t i = 0; i < (std::int64_t)n; ++i) {
    const double pa = (double)tot[(std::size_t)i] / U;
    for (std::size_t j = 0; j <= (std::size_t)i; ++j) {
      const double inter = (double)C[(std::size_t)i * n + j];
      const double uni = (double)tot[(std::size_t)i] + (double)tot[j] - inter;
      const double v = jaccard_value(inter, uni, pa, (double)tot[j] / U, centered);
      out[(std::size_t)i * n + j] = v;
      out[j * n + (std::size_t)i] = v;
    }
  }
}

} // namespace pam

// pam_binning.hpp -- index functions mapping environmental coordinates to bin
// ids. Two modes share the same downstream contract: a dense id in
// [0, n_bins()) or -1 for non-finite input.
//
//   GridBinning       -- uniform D-dimensional grid; the classical choice at
//                        D <= 3 and what c_bin_assign implements at D = 2.
//   DictionaryBinning -- realised environmental units: deduplicated quantised
//                        tuples get dense ids in insertion rank. This is the
//                        19-dimensional path: a grid nb^D is infeasible for
//                        D > 3, while the realised-unit count stays affordable
//                        for the PamBitset (and the bitset does not care how
//                        ids were produced).
#pragma once

#include <cstddef>
#include <cstdint>
#include <cmath>
#include <unordered_map>
#include <vector>

namespace pam {

class GridBinning {
public:
  // mins/maxs/nbs have length D. Bin id =
  //   sum_d clamp(floor((x_d - min_d) * inv_d), 0, nb_d - 1) * stride_d,
  // with stride_d = prod_{k > d} nb_k: dim 0 is the MOST significant digit,
  // matching c_bin_assign's id = bx * nb + by.
  GridBinning(std::vector<double> mins, std::vector<double> maxs,
              std::vector<int> nbs)
      : min_(std::move(mins)), nb_(std::move(nbs)) {
    const std::size_t D = min_.size();
    inv_.resize(D);
    stride_.resize(D);
    n_bins_ = 1;
    for (std::size_t d = 0; d < D; ++d) {
      inv_[d] = (double)nb_[d] / (maxs[d] - min_[d]);
      n_bins_ *= (std::size_t)nb_[d];
    }
    for (std::size_t d = 0; d < D; ++d) {
      std::size_t s = 1;
      for (std::size_t k = d + 1; k < D; ++k) s *= (std::size_t)nb_[k];
      stride_[d] = s;
    }
  }

  // x[d] read as x[d * step]: step = 1 for contiguous points, step = n_rows
  // for a row of a column-major matrix (zero-copy binning of R matrices).
  inline std::int64_t assign_strided(const double* x, std::ptrdiff_t step) const noexcept {
    std::int64_t id = 0;
    for (std::size_t d = 0; d < min_.size(); ++d) {
      const double xd = x[d * step];
      if (!std::isfinite(xd)) return -1;
      long v = (long)((xd - min_[d]) * inv_[d]);
      if (v < 0) v = 0; else if (v >= nb_[d]) v = nb_[d] - 1;
      id += (std::int64_t)v * (std::int64_t)stride_[d];
    }
    return id;
  }

  inline std::int64_t assign(const double* x) const noexcept {
    return assign_strided(x, 1);
  }

  std::size_t n_dims() const noexcept { return min_.size(); }
  std::size_t n_bins() const noexcept { return n_bins_; }

private:
  std::vector<double> min_, inv_;
  std::vector<int> nb_;
  std::vector<std::size_t> stride_;
  std::size_t n_bins_;
};

class DictionaryBinning {
public:
  // res[d] = quantisation width of dim d; tuples are floor(x_d / res_d).
  explicit DictionaryBinning(std::vector<double> resolutions)
      : res_(std::move(resolutions)), key_(res_.size()) {}

  // Insert-if-unseen dictionary assign. Ids are insertion-rank (order of
  // first sighting): dense and deterministic for a fixed input order, but
  // NOT stable across reorderings -- if stable ids are needed, deduplicate
  // in a sorted second pass instead.
  inline std::int64_t assign(const double* x) {
    for (std::size_t d = 0; d < res_.size(); ++d) {
      if (!std::isfinite(x[d])) return -1;
      key_[d] = (std::int32_t)std::floor(x[d] / res_[d]);
    }
    auto it = ids_.find(key_);
    if (it != ids_.end()) return it->second;
    const std::int64_t id = (std::int64_t)ids_.size();
    ids_.emplace(key_, id);
    return id;
  }

  std::size_t n_dims() const noexcept { return res_.size(); }
  std::size_t n_bins() const noexcept { return ids_.size(); }

private:
  struct KeyHash {
    std::size_t operator()(const std::vector<std::int32_t>& v) const noexcept {
      std::size_t h = v.size();
      for (std::int32_t k : v) h = h * 1315423911u ^ (std::size_t)(k + 0x9e3779b9u);
      return h;
    }
  };
  std::vector<double> res_;
  std::vector<std::int32_t> key_;   // per-instance scratch (single-threaded)
  std::unordered_map<std::vector<std::int32_t>, std::int64_t, KeyHash> ids_;
};

} // namespace pam

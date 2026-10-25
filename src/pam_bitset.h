// pam_bitset.h -- shared species x bin bitset accumulator.
//
// Pure C++17: no Rcpp/R headers, so the class is reusable by the engine
// (tile loop, rasteriser) without dragging Rcpp in, and is unit-testable
// outside R. One shared bitset for all threads; writes are lock-free atomic
// ORs, so memory is n_species*n_bins/8 bytes TOTAL and never scales with the
// thread count or the number of input pairs.
//
// Rows are padded to a whole number of 64-bit words (<= 63 wasted bits per
// species): row(s) then points at n_bins contiguous bits starting on a word
// boundary, which is what pairwise kernels (Jaccard AND/OR popcount) need.
//
// Portability note: __atomic_fetch_or and __builtin_ctzll are GCC/Clang
// builtins; Rtools on Windows ships GCC, so this matches the portability
// profile the kernel already had.
#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace pam {

class PamBitset {
public:
  PamBitset(std::size_t n_species, std::size_t n_bins)
      : n_species_(n_species), n_bins_(n_bins),
        wpr_((n_bins + 63u) / 64u),
        words_(n_species * wpr_, 0ULL) {}

  // Lock-free set bit (safe to call inside OpenMP regions).
  inline void add_atomic(std::size_t s, std::size_t b) noexcept {
    __atomic_fetch_or(&words_[s * wpr_ + (b >> 6)], 1ULL << (b & 63),
                      __ATOMIC_RELAXED);
  }

  // Serial set bit (single-threaded paths).
  inline void add(std::size_t s, std::size_t b) noexcept {
    words_[s * wpr_ + (b >> 6)] |= (1ULL << (b & 63));
  }

  // Word-aligned view of species row s: wpr_ words, low n_bins_ bits used
  // (tail bits of the last word are never set).
  inline const std::uint64_t* row(std::size_t s) const noexcept {
    return &words_[s * wpr_];
  }
  std::size_t words_per_row() const noexcept { return wpr_; }

  // Serial marginal pass; call only after all writes are complete (the end of
  // the parallel region's implicit barrier suffices -- atomic relaxed ORs are
  // all visible after it). Returns nnz; fills richness[b] and footprint[s].
  template <typename OutIt>
  std::size_t marginals(OutIt richness, OutIt footprint) const {
    std::size_t nnz = 0;
    for (std::size_t s = 0; s < n_species_; ++s) {
      const std::uint64_t* r = row(s);
      for (std::size_t w = 0; w < wpr_; ++w) {
        std::uint64_t word = r[w];
        while (word) {
          const std::size_t b = w * 64u + (std::size_t)__builtin_ctzll(word);
          ++richness[b];
          ++footprint[s];
          ++nnz;
          word &= word - 1;
        }
      }
    }
    return nnz;
  }

  std::size_t n_species() const noexcept { return n_species_; }
  std::size_t n_bins() const noexcept { return n_bins_; }
  std::size_t bytes() const noexcept { return words_.size() * sizeof(std::uint64_t); }

  static std::size_t bytes_for(std::size_t n_species, std::size_t n_bins) noexcept {
    return n_species * ((n_bins + 63u) / 64u) * sizeof(std::uint64_t);
  }

private:
  std::size_t n_species_, n_bins_, wpr_;
  std::vector<std::uint64_t> words_;
};

} // namespace pam

// PamEngine.cpp -- tiled environmental-PAM engine.
//
// Pipeline per strip (full-width rows, so raster cells form a contiguous
// range):
//   EnvReader.read_strip -> GridBinning.assign per pixel -> availability
//   counts -> scatter of the strip's (species, cell) pairs into PamBitset
//   (lock-free atomic ORs) or a sort-fallback key vector.
//
// Pairs are sorted by cell once, so the pair stream and the strip stream
// march in a single merge pass -- no per-tile lookups and no dense
// intermediate of length P anywhere.
//
// The XPtr wrappers at the bottom are always compiled; GDAL absence surfaces
// only as a runtime error from compute().

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "PamEngine.h"
#include "xpam_gdal_io.h"
#include "pam_bitset.h"
#include "pam_binning.h"
#include "pam_jaccard.h"
#include "pam_occ_io.h"
#include "pam_rasterise.h"
#include "pam_threads.h"

using namespace Rcpp;

namespace pam {

void PamEngine::open_env(const std::vector<std::string>& files) {
  env_files_ = files;
}

void PamEngine::set_occ(Rcpp::IntegerVector species, Rcpp::IntegerVector cell,
                        int n_species) {
  if (species.size() != cell.size())
    throw std::runtime_error("PamEngine::set_occ: species and cell must have the same length");
  occ_sp_ = species;
  occ_cell_ = cell;
  n_species_ = n_species;
}

void PamEngine::set_bins(std::vector<double> mins, std::vector<double> maxs,
                         std::vector<int> nbs) {
  if (mins.size() != maxs.size() || mins.size() != nbs.size())
    throw std::runtime_error("PamEngine::set_bins: mins, maxs and nbs must have equal length");
  mins_ = std::move(mins); maxs_ = std::move(maxs); nbs_ = std::move(nbs);
}

void PamEngine::set_occ_file(const std::string& path) { occ_file_ = path; }

void PamEngine::set_occ_polygons(const std::vector<std::string>& paths,
                                 const std::vector<std::string>& layers,
                                 const std::string& field,
                                 const std::vector<std::string>& levels,
                                 bool all_touched, int chunk_rows,
                                 const std::string& dump_path) {
  poly_paths_ = paths;
  poly_layers_ = layers;
  poly_field_ = field;
  poly_levels_ = levels;
  poly_all_touched_ = all_touched;
  poly_chunk_rows_ = chunk_rows < 1 ? 256 : chunk_rows;
  poly_dump_ = dump_path;
}

void PamEngine::set_threads(int n)      { n_threads_ = pam::safe_threads(n); }
void PamEngine::set_tile_rows(int rows) { tile_rows_ = rows < 1 ? 1 : rows; }
void PamEngine::set_budget_mb(double mb){ budget_mb_ = mb < 0 ? 0 : mb; }
void PamEngine::set_space(int s)        { space_ = s == 1 ? 1 : 0; }
void PamEngine::set_jaccard(int mode)   { jaccard_ = (mode >= 0 && mode <= 2) ? mode : 0; }
void PamEngine::set_by(int b)           { by_ = b == 1 ? 1 : 0; }
void PamEngine::set_out_budget_mb(double mb) { out_budget_mb_ = mb < 0 ? 0 : mb; }
void PamEngine::set_jaccard_threshold(double tau) {
  threshold_ = (tau > 0.0 && tau <= 1.0) ? tau : 0.0;
}
void PamEngine::set_export_incidence(bool v) { export_incidence_ = v; }

void PamEngine::compute() {
  if (env_files_.empty())
    throw std::runtime_error("PamEngine: no env files configured (open_env first)");
  const bool geo = (space_ == 1);
  if (!geo && mins_.empty())
    throw std::runtime_error("PamEngine: no binning configured (set_bins first)");
  const bool have_mem  = occ_sp_.size() > 0;
  const bool have_file = !occ_file_.empty();
  const bool have_poly = !poly_paths_.empty();
  if (!have_mem && !have_file && !have_poly)
    throw std::runtime_error("PamEngine: no occurrence source configured");
  if (have_mem && n_species_ <= 0)
    throw std::runtime_error("PamEngine: n_species must be positive");

  EnvReader env(env_files_);          // throws when built without GDAL
  const int nc = env.ncols(), nr = env.nrows(), D = env.n_dims();
  if (!geo && (int)mins_.size() != D)
    throw std::runtime_error("PamEngine: bins dimension != env raster dims");

  GridBinning gb(mins_, maxs_, nbs_);   // unused when geo (bin id = cell id)
  const bool by_site = (by_ == 1);
  const std::int64_t ncells = (std::int64_t)nc * nr;
  const std::int64_t n_bins = geo ? ncells : (std::int64_t)gb.n_bins();
  res_.ncols = nc; res_.nrows = nr;
  res_.space = geo ? 1 : 0;
  res_.jaccard_mode = jaccard_;
  res_.by = by_site ? 1 : 0;
  res_.threshold = threshold_;
  res_.export_incidence = export_incidence_;

  // Occurrence source: streamed file (Phase C) or SEXP-anchored memory.
  // Both yield (cell, sp) pairs sorted by cell exactly once; the file path
  // additionally bounds its read window, so the pair table never
  // materializes in either language.
  std::int64_t oor = 0;
  std::unique_ptr<OccSource> src;
  if (have_file) {
    src.reset(new FileOccSource(occ_file_, 1u << 22, &oor, ncells));
    n_species_ = (int)src->n_species();
  } else if (have_poly) {
    // The polygon grid IS the env raster's: geotransform + SRS come from
    // the same EnvReader that feeds the strips.
    double gt[6];
    env.geotransform(gt);
    if (poly_paths_.size() > 1 && poly_levels_.empty())
      throw std::runtime_error(
          "PamEngine: multiple polygon datasets need explicit `species` "
          "levels (species codes must be consistent across files)");
    std::vector<std::unique_ptr<OccSource>> polys;
    polys.reserve(poly_paths_.size());
    for (std::size_t k = 0; k < poly_paths_.size(); ++k) {
      const std::string layer =
          k < poly_layers_.size() ? poly_layers_[k] : std::string();
      polys.emplace_back(new PolygonOccSource(
          poly_paths_[k], layer, poly_field_, poly_levels_,
          nc, nr, gt, env.srs_wkt(), poly_all_touched_, poly_chunk_rows_));
    }
    n_species_ = (int)polys[0]->n_species();
    if (n_species_ <= 0)
      throw std::runtime_error(
          "PamEngine: polygon source produced no species");
    auto* first = static_cast<PolygonOccSource*>(polys[0].get());
    res_.species = first->species_levels();
    src = polys.size() == 1
        ? std::move(polys[0])
        : std::unique_ptr<OccSource>(new MergeOccSource(std::move(polys)));
    // The dump tees the stream the engine actually consumes -- merged
    // and deduplicated across datasets.
    if (!poly_dump_.empty())
      src.reset(new DumpOccSource(std::move(src), poly_dump_));
  } else {
    src.reset(new MemOccSource(INTEGER(occ_sp_), INTEGER(occ_cell_),
                               occ_sp_.size(), n_species_, ncells, &oor));
  }

  // Accumulator geometry: species mode packs rows = species, bits = units;
  // site mode transposes (rows = units, bits = species) so the same row
  // kernels yield site x site Jaccard.  Site mode has no fallbacks: the
  // co-occurrence path would need a species-sorted stream and the sort
  // kernel keys on species -- the transposed bitset must fit the budget.
  const std::size_t bs_rows = by_site ? (std::size_t)n_bins
                                      : (std::size_t)n_species_;
  const std::size_t bs_cols = by_site ? (std::size_t)n_species_
                                      : (std::size_t)n_bins;
  res_.state_mb = (double)PamBitset::bytes_for(bs_rows, bs_cols) / 1048576.0;
  const bool use_bitset = res_.state_mb <= budget_mb_;
  if (by_site && !use_bitset)
    throw std::runtime_error(
        "PamEngine: by='site' needs the transposed bitset (units x species, "
        + std::to_string(res_.state_mb) + " MB) within max_bitset_mb; no "
        "co-occurrence fallback exists (it requires a species-sorted "
        "stream). Reduce the grid/bins or raise max_bitset_mb");
  const bool use_cooc = geo && !by_site && !use_bitset;
  res_.kernel = use_bitset ? "bitmap" : (geo ? "cooc" : "sort");
  if (export_incidence_ && use_cooc)
    throw std::runtime_error(
        "PamEngine: incidence export needs the bitset or sort accumulator; "
        "the geo co-occurrence path does not retain (species, cell) pairs. "
        "Raise max_bitset_mb, reduce the grid, or use occ_dump (XPOC) for "
        "the on-disk sparse form");

  // The Jaccard result itself is a dense n^2 double matrix; refuse early
  // when it cannot be marshalled to R (e.g. site mode over a large grid).
  // Bounded by out_budget_mb_, not max_bitset_mb: one governs the
  // incidence state, the other the result handed back to R.
  if (jaccard_ != 0 && threshold_ <= 0.0) {
    const std::size_t n_ent = by_site ? (std::size_t)n_bins
                                      : (std::size_t)n_species_;
    const double out_mb = (double)n_ent * (double)n_ent * 8.0 / 1048576.0;
    if (out_mb > out_budget_mb_)
      throw std::runtime_error(
          "PamEngine: dense Jaccard output would be "
          + std::to_string(out_mb) + " MB > max_out_mb. For "
          + (by_site ? "site mode" : "species mode")
          + " over this many entities use a smaller grid/binning "
          + "(a thresholded edge list is planned but not implemented)");
  }

  std::unique_ptr<PamBitset> bs;
  if (use_bitset)
    bs.reset(new PamBitset(bs_rows, bs_cols));
  std::vector<std::uint64_t> acc;      // sort-fallback: (sp << 32) | bin
  const std::size_t nsp = (std::size_t)n_species_;
  if (use_cooc)
    res_.cooc.assign(nsp * nsp, 0);    // C[a*n+b] = cells shared by a and b

  res_.richness.assign(use_cooc ? 0 : (std::size_t)n_bins, 0);
  res_.footprint.assign(nsp, 0);
  res_.availability.assign(geo ? 0 : (std::size_t)n_bins, 0.0);

  const int nt = n_threads_ > 0 ? n_threads_ : 1;
  std::vector<double> block;
  std::vector<std::int32_t> binbuf;

  // Single merge pass: the source yields pairs sorted by cell while strips
  // advance in cell order, so each pair is seen exactly once.  Per strip we
  // drain the pairs it owns into a small buffer, then scatter in parallel.
  std::vector<std::int64_t> pc;
  std::vector<std::int32_t> ps;
  std::int64_t cur_cell;  std::int32_t cur_sp;
  bool have = src->next(cur_cell, cur_sp);

  for (int y0 = 0; y0 < nr; y0 += tile_rows_) {
    const int h = std::min(tile_rows_, nr - y0);
    const std::int64_t npix = (std::int64_t)h * nc;
    const std::int64_t lo = (std::int64_t)y0 * nc;
    const std::int64_t hi = lo + npix;

    if (!geo) {
      block.resize((std::size_t)npix * D);
      binbuf.resize((std::size_t)npix);
      env.read_strip(y0, h, block.data());                     // serial I/O

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nt)
#endif
      for (std::int64_t i = 0; i < npix; ++i)
        binbuf[(std::size_t)i] = (std::int32_t)gb.assign(&block[(std::size_t)i * D]);

      for (std::int64_t i = 0; i < npix; ++i) {
        const std::int32_t b = binbuf[(std::size_t)i];
        if (b >= 0) res_.availability[(std::size_t)b] += 1.0;
      }
    }

    pc.clear(); ps.clear();
    while (have && cur_cell < hi) {          // pairs are sorted: < lo means
      if (cur_cell >= lo) {                  // they belong to a past strip --
        pc.push_back(cur_cell);              // cannot happen post-filter, but
        ps.push_back(cur_sp);                // skip defensively below
      }
      have = src->next(cur_cell, cur_sp);
    }

    if (use_bitset) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nt)
#endif
      for (std::int64_t k = 0; k < (std::int64_t)pc.size(); ++k) {
        const std::int64_t b = geo ? pc[(std::size_t)k]
          : (std::int64_t)binbuf[(std::size_t)(pc[(std::size_t)k] - lo)];
        if (b >= 0 && b < n_bins) {
          if (by_site)
            bs->add_atomic((std::size_t)b, (std::size_t)ps[(std::size_t)k]);
          else
            bs->add_atomic((std::size_t)ps[(std::size_t)k], (std::size_t)b);
        }
      }
    } else if (use_cooc) {
      // Cell-sorted runs: each run is the species set of one cell.
      // Dedup within the run (mem sources may carry duplicate pairs),
      // then accumulate unordered pairs into the co-occurrence matrix
      // and per-species cell counts.  Runs are independent; updates are
      // atomic adds over a shared n^2 table.
      std::vector<std::int64_t> starts;
      for (std::size_t k = 0; k < pc.size();) {
        starts.push_back((std::int64_t)k);
        std::size_t k1 = k + 1;
        while (k1 < pc.size() && pc[k1] == pc[k]) ++k1;
        k = k1;
      }
#ifdef _OPENMP
#pragma omp parallel num_threads(nt)
#endif
      {
        std::vector<std::int32_t> run;
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 16)
#endif
        for (std::size_t r = 0; r < starts.size(); ++r) {
          const std::size_t k0 = (std::size_t)starts[r];
          const std::size_t k1 = r + 1 < starts.size()
            ? (std::size_t)starts[r + 1] : pc.size();
          run.assign(ps.begin() + (std::ptrdiff_t)k0,
                     ps.begin() + (std::ptrdiff_t)k1);
          std::sort(run.begin(), run.end());
          run.erase(std::unique(run.begin(), run.end()), run.end());
          for (std::size_t a = 0; a < run.size(); ++a) {
            const std::size_t sa = (std::size_t)run[a];
            __atomic_fetch_add(&res_.footprint[sa], 1, __ATOMIC_RELAXED);
            for (std::size_t bb = a + 1; bb < run.size(); ++bb) {
              const std::size_t sb = (std::size_t)run[bb];
              __atomic_fetch_add(&res_.cooc[sa * nsp + sb], 1,
                                 __ATOMIC_RELAXED);
              __atomic_fetch_add(&res_.cooc[sb * nsp + sa], 1,
                                 __ATOMIC_RELAXED);
            }
          }
        }
      }
    } else {
      for (std::size_t k = 0; k < pc.size(); ++k) {
        const std::int32_t b = binbuf[(std::size_t)(pc[k] - lo)];
        if (b >= 0)
          acc.push_back(((std::uint64_t)(std::uint32_t)ps[k] << 32)
                        | (std::uint32_t)b);
      }
    }
    res_.n_strips++;
  }

  if (use_bitset) {
    // marginals(arg1, arg2): arg1 <- per-column counts, arg2 <- per-row
    // counts.  In the transposed bitset the roles swap, so pass the
    // result vectors in the order that keeps the public semantics:
    // richness = species per unit, footprint = units per species.
    if (by_site)
      res_.nnz = (std::int64_t)bs->marginals(res_.footprint.data(),
                                             res_.richness.data());
    else
      res_.nnz = (std::int64_t)bs->marginals(res_.richness.data(),
                                             res_.footprint.data());
  } else if (use_cooc) {
    for (std::size_t s = 0; s < nsp; ++s) {
      res_.cooc[s * nsp + s] = (std::uint64_t)res_.footprint[s]; // diagonal
      res_.nnz += res_.footprint[s];
    }
  } else {
    std::sort(acc.begin(), acc.end());
    std::int64_t nnz = 0;
    for (std::size_t k = 0; k < acc.size(); ++k) {
      if (k > 0 && acc[k] == acc[k - 1]) continue;
      const std::int64_t s = (std::int64_t)(acc[k] >> 32);
      const std::int64_t b = (std::int64_t)(acc[k] & 0xffffffffu);
      res_.richness[(std::size_t)b]++; res_.footprint[(std::size_t)s]++; nnz++;
    }
    res_.nnz = nnz;
  }
  res_.n_out_of_range = oor;
  res_.n_skipped = src->n_skipped();

  // Optional incidence export: sparse triplets (i = unit, j = species).
  // Two int32 vectors of nnz elements must fit out_budget_mb_; beyond
  // that the honest shape is the XPOC occ_dump on disk, not an R object.
  if (export_incidence_) {
    const double inc_mb = (double)res_.nnz * 16.0 / 1048576.0;
    if (inc_mb > out_budget_mb_)
      throw std::runtime_error(
          "PamEngine: incidence has nnz=" + std::to_string(res_.nnz) +
          " (~" + std::to_string(inc_mb) +
          " MB as triplets) > max_out_mb; write it to disk via occ_dump "
          "instead of returning it to R");
    res_.inc_i.reserve((std::size_t)res_.nnz);
    res_.inc_j.reserve((std::size_t)res_.nnz);
    if (use_bitset) {
      for (std::size_t s = 0; s < bs->n_species(); ++s) {
        const std::uint64_t* r = bs->row(s);
        for (std::size_t w = 0; w < bs->words_per_row(); ++w) {
          std::uint64_t word = r[w];
          while (word) {
            const std::int32_t b =
                (std::int32_t)(w * 64u + (std::size_t)__builtin_ctzll(word));
            if (by_site) {
              res_.inc_i.push_back((std::int32_t)s);
              res_.inc_j.push_back(b);
            } else {
              res_.inc_i.push_back(b);
              res_.inc_j.push_back((std::int32_t)s);
            }
            word &= word - 1;
          }
        }
      }
    } else {
      std::int64_t prev = -1;
      for (std::size_t k = 0; k < acc.size(); ++k) {
        if ((std::int64_t)acc[k] == prev) continue;   // dedup same as nnz pass
        prev = (std::int64_t)acc[k];
        res_.inc_i.push_back((std::int32_t)(acc[k] & 0xffffffffu));
        res_.inc_j.push_back((std::int32_t)(acc[k] >> 32));
      }
    }
  }

  // Optional Jaccard post-pass over whichever accumulator ran
  // (fastJaccard semantics: J = |n|/|u|, expected value on empty unions,
  // optional centring).  With threshold_ > 0 the output is the exact
  // edge list via a prefix-filtered join instead of the dense matrix.
  if (jaccard_ != 0) {
    const bool centered = (jaccard_ == 2);
    if (threshold_ > 0.0) {
      if (use_cooc) {
        pam::jaccard_edges_from_cooc(res_.cooc, res_.footprint, nsp,
                                     (std::size_t)n_bins, threshold_,
                                     centered, nt,
                                     res_.edge_i, res_.edge_j, res_.edge_v);
      } else {
        // Unit-id lists per entity: from the bitset rows (either
        // orientation -- rows are the entities) or the sort-fallback
        // accumulator.  freq[] for the prefix choice is the global
        // per-column count: richness in species mode, footprint in
        // site mode.
        std::vector<std::vector<std::int32_t>> sets(
            use_bitset ? bs->n_species() : nsp);
        if (use_bitset) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nt)
#endif
          for (std::int64_t s = 0; s < (std::int64_t)bs->n_species(); ++s) {
            const std::uint64_t* r = bs->row((std::size_t)s);
            std::vector<std::int32_t>& v = sets[(std::size_t)s];
            for (std::size_t w = 0; w < bs->words_per_row(); ++w) {
              std::uint64_t word = r[w];
              while (word) {
                v.push_back((std::int32_t)(w * 64u +
                    (std::size_t)__builtin_ctzll(word)));
                word &= word - 1;
              }
            }
          }
        } else {
          for (std::size_t k = 0; k < acc.size(); ++k) {
            std::vector<std::int32_t>& s = sets[(std::size_t)(acc[k] >> 32)];
            const std::int32_t b = (std::int32_t)(acc[k] & 0xffffffffu);
            if (s.empty() || s.back() != b) s.push_back(b);
          }
        }
        const std::vector<std::int64_t>& freq =
            by_site ? res_.footprint : res_.richness;
        const std::size_t U =
            by_site ? nsp : (std::size_t)n_bins;
        pam::jaccard_edges(sets, freq, U, threshold_, centered, nt,
                           res_.edge_i, res_.edge_j, res_.edge_v);
      }
    } else if (use_bitset) {
      pam::jaccard_rows(*bs, centered, nt, res_.jaccard);
    } else if (use_cooc) {
      pam::jaccard_from_cooc(res_.cooc, res_.footprint, nsp,
                             (std::size_t)n_bins, centered, nt,
                             res_.jaccard);
    } else {
      // acc may still carry duplicate (sp, bin) pairs (the nnz pass above
      // skips them without erasing); dedup while building the lists or the
      // union count is inflated.
      std::vector<std::vector<std::int32_t>> sets(nsp);
      for (std::size_t k = 0; k < acc.size(); ++k) {
        std::vector<std::int32_t>& s = sets[(std::size_t)(acc[k] >> 32)];
        const std::int32_t b = (std::int32_t)(acc[k] & 0xffffffffu);
        if (s.empty() || s.back() != b) s.push_back(b);
      }
      pam::jaccard_lists(sets, (std::size_t)n_bins, centered, nt,
                         res_.jaccard);
    }
  }
}

} // namespace pam

// ── XPtr wrappers (always compiled) ─────────────────────────────────────────

//' Create a new PamEngine instance
//' @return An \code{externalptr} to a new C++ \code{PamEngine} object.
//' @noRd
//' @noRd
// [[Rcpp::export]]
SEXP xpam_engine_create() {
  Rcpp::XPtr<pam::PamEngine> ptr(new pam::PamEngine(), true);
  return ptr;
}

//' Set env raster files (1 multi-band or D single-band)
//' @noRd
// [[Rcpp::export]]
void xpam_engine_open_env(SEXP xptr, std::vector<std::string> files) {
  Rcpp::XPtr<pam::PamEngine>(xptr)->open_env(files);
}

//' Attach the (species, cell) pair table (0-based codes, zero-copy)
//' @noRd
// [[Rcpp::export]]
void xpam_engine_set_occ(SEXP xptr, Rcpp::IntegerVector species,
                         Rcpp::IntegerVector cell, int n_species) {
  Rcpp::XPtr<pam::PamEngine>(xptr)->set_occ(species, cell, n_species);
}

//' Attach a streamed occurrence file (XPOC v1, cell-sorted)
//' @noRd
// [[Rcpp::export]]
void xpam_engine_set_occ_file(SEXP xptr, std::string path) {
  Rcpp::XPtr<pam::PamEngine>(xptr)->set_occ_file(path);
}

//' Attach polygon occurrence source(s): parallel vectors of dataset
//' paths and layer names ("" = first layer), shared species field, code
//' levels (empty = derive; required when paths > 1), ALL_TOUCHED flag,
//' rasterise chunk height, and an optional XPOC dump path
//' @noRd
// [[Rcpp::export]]
void xpam_engine_set_occ_polygons(SEXP xptr, std::vector<std::string> paths,
                                  std::vector<std::string> layers,
                                  std::string field,
                                  Rcpp::CharacterVector levels,
                                  bool all_touched, int chunk_rows,
                                  std::string dump_path) {
  Rcpp::XPtr<pam::PamEngine>(xptr)->set_occ_polygons(
      paths, layers, field, Rcpp::as<std::vector<std::string>>(levels),
      all_touched, chunk_rows, dump_path);
}

//' Set the D-dimensional grid binning (mins, maxs, nbs of length D)
//' @noRd
// [[Rcpp::export]]
void xpam_engine_set_bins(SEXP xptr, Rcpp::NumericVector mins,
                          Rcpp::NumericVector maxs, Rcpp::IntegerVector nbs) {
  Rcpp::XPtr<pam::PamEngine>(xptr)->set_bins(
    Rcpp::as<std::vector<double>>(mins), Rcpp::as<std::vector<double>>(maxs),
    Rcpp::as<std::vector<int>>(nbs));
}

//' Set OpenMP threads (capped by OMP_THREAD_LIMIT)
//' @noRd
// [[Rcpp::export]]
void xpam_engine_set_threads(SEXP xptr, int n) {
  Rcpp::XPtr<pam::PamEngine>(xptr)->set_threads(n);
}

//' Set strip height in rows (default 512)
//' @noRd
// [[Rcpp::export]]
void xpam_engine_set_tile_rows(SEXP xptr, int rows) {
  Rcpp::XPtr<pam::PamEngine>(xptr)->set_tile_rows(rows);
}

//' Set the bitset memory budget in MiB (sort fallback beyond it)
//' @noRd
// [[Rcpp::export]]
void xpam_engine_set_budget_mb(SEXP xptr, double mb) {
  Rcpp::XPtr<pam::PamEngine>(xptr)->set_budget_mb(mb);
}

//' Set the incidence space: 0 = env bins, 1 = geographic cells
//' @noRd
// [[Rcpp::export]]
void xpam_engine_set_space(SEXP xptr, int space) {
  Rcpp::XPtr<pam::PamEngine>(xptr)->set_space(space);
}

//' Set the Jaccard post-pass: 0 = none, 1 = raw, 2 = centered
//' @noRd
// [[Rcpp::export]]
void xpam_engine_set_jaccard(SEXP xptr, int mode) {
  Rcpp::XPtr<pam::PamEngine>(xptr)->set_jaccard(mode);
}

//' Set the Jaccard entities: 0 = species, 1 = site (transposed)
//' @noRd
// [[Rcpp::export]]
void xpam_engine_set_by(SEXP xptr, int by) {
  Rcpp::XPtr<pam::PamEngine>(xptr)->set_by(by);
}

//' Set the dense-output feasibility cap for the Jaccard matrix (MB)
//' @noRd
// [[Rcpp::export]]
void xpam_engine_set_out_budget_mb(SEXP xptr, double mb) {
  Rcpp::XPtr<pam::PamEngine>(xptr)->set_out_budget_mb(mb);
}

//' Set the Jaccard edge threshold in (0,1]; 0 = dense matrix output
//' @noRd
// [[Rcpp::export]]
void xpam_engine_set_jaccard_threshold(SEXP xptr, double tau) {
  Rcpp::XPtr<pam::PamEngine>(xptr)->set_jaccard_threshold(tau);
}

//' Export the incidence as sparse triplets (i = unit, j = species)
//' @noRd
// [[Rcpp::export]]
void xpam_engine_set_export_incidence(SEXP xptr, bool v) {
  Rcpp::XPtr<pam::PamEngine>(xptr)->set_export_incidence(v);
}

//' Run the engine and return the result list
//' @noRd
// [[Rcpp::export]]
Rcpp::List xpam_engine_compute(SEXP xptr) {
  Rcpp::XPtr<pam::PamEngine> eng(xptr);
  eng->compute();
  const pam::PamResult& r = eng->results();
  Rcpp::IntegerVector rich(r.richness.size()), foot(r.footprint.size());
  for (std::size_t i = 0; i < r.richness.size(); ++i) rich[i] = (int)r.richness[i];
  for (std::size_t i = 0; i < r.footprint.size(); ++i) foot[i] = (int)r.footprint[i];
  Rcpp::List out = Rcpp::List::create(
    _["nnz"] = (double)r.nnz,
    _["richness"] = rich,
    _["footprint"] = foot,
    _["availability"] = Rcpp::wrap(r.availability),
    _["kernel"] = r.kernel,
    _["state_mb"] = r.state_mb,
    _["n_out_of_range"] = (double)r.n_out_of_range,
    _["n_skipped"] = (double)r.n_skipped,
    _["species_levels"] = Rcpp::wrap(r.species),
    _["ncols"] = r.ncols, _["nrows"] = r.nrows, _["n_strips"] = r.n_strips,
    _["space"] = r.space, _["jaccard_mode"] = r.jaccard_mode,
    _["by"] = r.by, _["threshold"] = r.threshold);
  if (!r.jaccard.empty()) {
    const int n = (int)std::lround(std::sqrt((double)r.jaccard.size()));
    Rcpp::NumericMatrix jm(n, n);
    std::copy(r.jaccard.begin(), r.jaccard.end(), jm.begin());
    out["jaccard"] = jm;
  }
  if (!r.cooc.empty()) {
    const int n = (int)std::lround(std::sqrt((double)r.cooc.size()));
    Rcpp::NumericMatrix cm(n, n);
    for (std::size_t k = 0; k < r.cooc.size(); ++k)
      cm[k] = (double)r.cooc[k];
    out["cooc"] = cm;
  }
  if (r.threshold > 0.0) {           // edge mode: marshal even when empty
    const std::size_t ne = r.edge_i.size();
    Rcpp::IntegerVector ii(ne), jj(ne);
    for (std::size_t k = 0; k < ne; ++k) {   // 0-based -> 1-based for R
      ii[k] = r.edge_i[k] + 1;
      jj[k] = r.edge_j[k] + 1;
    }
    out["edge_i"] = ii;
    out["edge_j"] = jj;
    out["edge_v"] = Rcpp::NumericVector(r.edge_v.begin(), r.edge_v.end());
  }
  if (r.export_incidence) {            // marshal even when empty
    const std::size_t nz = r.inc_i.size();
    Rcpp::IntegerVector ii(nz), jj(nz);
    for (std::size_t k = 0; k < nz; ++k) {   // 0-based -> 1-based for R
      ii[k] = r.inc_i[k] + 1;
      jj[k] = r.inc_j[k] + 1;
    }
    out["inc_i"] = ii;
    out["inc_j"] = jj;
  }
  return out;
}

//' Was the package built with GDAL? (compile-time flag)
//' @noRd
// [[Rcpp::export]]
bool xpam_have_gdal() {
  return pam::EnvReader::built_with_gdal();
}

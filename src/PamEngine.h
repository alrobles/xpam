// PamEngine.h -- tiled environmental-PAM engine for xpam.
//
// Orchestrates the full pipeline, mirroring xbioclim::BioclimEngine:
//   open_env() -> set_occ() -> set_bins() -> [set_threads() / set_tile_rows()
//   / set_budget_mb()] -> compute()
//
// All raster I/O and computation live in C++; R sends only file paths, the
// (species, cell) pair table, and scalar parameters.  Peak memory is
// O(bitset + strip), never O(n_cells) env layers in R and never a dense
// species-by-site matrix.
//
// Cell convention: occ cells are 0-based raster cell indexes
// (row * ncols + col), matching GDAL/terra cellnumber minus 1.  The R
// conductor converts terra cellnumbers before calling set_occ().
//
// compute() requires GDAL; without it every EnvReader call throws a clear
// std::runtime_error (re-raised as Rcpp::stop by the XPtr wrapper).
#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>
#include <Rcpp.h>

namespace pam {

struct PamResult {
  std::vector<int64_t> richness;       // species per bin (env; geo: per cell)
  std::vector<int64_t> footprint;      // bins per species (geo: cells)
  std::vector<double>  availability;   // cells per bin (env space only)
  std::vector<std::string> species;    // levels of a polygon source (else empty)
  std::vector<double>  jaccard;        // n x n symmetric similarity (on request)
  std::vector<std::int32_t> edge_i, edge_j;   // thresholded edges (i<j)
  std::vector<double>  edge_v;         //   edge value (raw or centred J)
  std::vector<std::int32_t> inc_i, inc_j;     // incidence triplets (unit, sp)
  double threshold = 0.0;              // active Jaccard threshold (0 = dense)
  bool export_incidence = false;       // incidence export was requested
  std::vector<std::uint64_t> cooc;     // n x n co-occurrence counts (geo fallback)
  std::int64_t nnz = 0;
  std::int64_t n_out_of_range = 0;
  std::int64_t n_skipped = 0;          // inputs unusable (e.g. bad field/geom)
  std::string kernel;
  int space = 0;                       // 0 = env bins, 1 = geographic cells
  int jaccard_mode = 0;                // 0 none, 1 raw, 2 centered
  int by = 0;                          // 0 = species, 1 = site (transposed)
  double state_mb = 0.0;
  int ncols = 0, nrows = 0, n_strips = 0;
};

class PamEngine {
public:
  PamEngine()  = default;
  ~PamEngine() = default;

  PamEngine(const PamEngine&)            = delete;
  PamEngine& operator=(const PamEngine&) = delete;
  PamEngine(PamEngine&&)                 = delete;
  PamEngine& operator=(PamEngine&&)      = delete;

  // env raster paths: one multi-band file (D = nbands) or D single-band files.
  void open_env(const std::vector<std::string>& files);

  // (species, cell) pair table; species/cell are 0-based codes.
  // The vectors are stored by reference (SEXP-anchored, zero-copy).
  void set_occ(Rcpp::IntegerVector species, Rcpp::IntegerVector cell,
               int n_species);

  // Occurrence file (XPOC v1, cell-sorted); streamed in chunks instead of
  // materializing the pair table.  Overrides set_occ when non-empty.
  void set_occ_file(const std::string& path);

  // Polygon occurrence source (Phase D): rasterise the `field` attribute of
  // a vector layer onto the env grid, streamed chunk-wise through the same
  // merge pass.  `levels` orders species codes (empty = file order);
  // `dump_path` optionally writes the emitted pairs as an XPOC file.
  // Precedence over set_occ/set_occ_file matches set_occ_file.
  // Polygon source(s).  `paths` and `layers` are parallel vectors
  // (layers may be empty -> first layer of each dataset); more than one
  // dataset is k-way merged by cell, and then `levels` must be given
  // (species codes have to be consistent across datasets).
  void set_occ_polygons(const std::vector<std::string>& paths,
                        const std::vector<std::string>& layers,
                        const std::string& field,
                        const std::vector<std::string>& levels,
                        bool all_touched, int chunk_rows,
                        const std::string& dump_path);

  // Grid binning in D dims (length must equal env dims at compute()).
  void set_bins(std::vector<double> mins, std::vector<double> maxs,
                std::vector<int> nbs);

  void set_threads(int n);
  void set_tile_rows(int rows);      // rows per strip; default 512
  void set_budget_mb(double mb);     // bitset budget; default 1024

  // Incidence unit space: 0 = environmental bins (GridBinning over env
  // values), 1 = geographic cells (bin id == cell index; env files then
  // only provide grid geometry -- values are never read).  When the
  // species x cell bitset would exceed the budget, geo space switches to
  // an O(n_species^2) co-occurrence accumulator fed by the same
  // cell-sorted pair stream.
  void set_space(int s);
  // Pairwise Jaccard post-pass over the accumulated incidence:
  // 0 = none, 1 = raw, 2 = centered (subtracts the probabilistic
  // expectation pa*pb/(pa+pb-pa*pb), fastJaccard semantics).
  void set_jaccard(int mode);
  // Jaccard entities: 0 = species (rows = species, bits = units; default),
  // 1 = site (transposed: rows = units, bits = species -> J between units
  // given the species they share).  by = site requires the transposed
  // bitset to fit the budget -- there is no co-occurrence fallback (that
  // path would need a species-sorted pair stream) and no sort fallback.
  void set_by(int b);
  // Feasibility cap for the dense n^2 Jaccard return (MB; default 4096).
  // Separate from max_bitset_mb: the bitset budget bounds the incidence
  // state, this bounds the result matrix that gets marshalled to R.
  void set_out_budget_mb(double mb);
  // Jaccard output threshold in (0, 1]: when > 0 the post-pass runs an
  // exact prefix-filtered join (AllPairs-style) and returns the edge
  // list (i, j, J >= tau) instead of the dense matrix -- output-
  // sensitive work, and the max_out_mb dense cap does not apply.
  void set_jaccard_threshold(double tau);
  // When true, export the incidence itself as sparse triplets
  // (i = unit, j = species; 0-based internally).  Requires the bitset or
  // sort accumulator -- the geo co-occurrence path does not retain
  // incidence -- and nnz * 16 B <= out_budget_mb_ to marshal to R.  For
  // larger streams the on-disk answer is the XPOC occ_dump, not an R
  // object.
  void set_export_incidence(bool v);

  // Run the pipeline and fill results().  Throws std::runtime_error when:
  //   * GDAL is not compiled in            * env dims != bins dims
  //   * inputs not configured              * any file cannot be read
  void compute();

  const PamResult& results() const { return res_; }

private:
  std::vector<std::string> env_files_;
  Rcpp::IntegerVector occ_sp_, occ_cell_;   // SEXP-anchored, zero-copy
  std::string occ_file_;                    // streamed source (Phase C)
  std::vector<std::string> poly_paths_, poly_layers_;              // Phase D+
  std::string poly_field_, poly_dump_;
  std::vector<std::string> poly_levels_;
  bool poly_all_touched_ = true;
  int poly_chunk_rows_ = 256;
  int n_species_ = 0;
  std::vector<double> mins_, maxs_;
  std::vector<int> nbs_;
  int n_threads_ = 1;
  int tile_rows_ = 512;
  double budget_mb_ = 1024.0;
  int space_ = 0;
  int jaccard_ = 0;
  int by_ = 0;
  double out_budget_mb_ = 4096.0;
  double threshold_ = 0.0;
  bool export_incidence_ = false;
  PamResult res_;
};

} // namespace pam

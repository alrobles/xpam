# xpam 1.0.0

First CRAN release.  `xpam` consolidates the GitHub-only 0.1.0--0.3.0
line into a single, reviewed submission: the small-input dense path
(`get_geo_pam()`) plus the streaming native engine, the sparse incidence
surface and the native Jaccard backend described under 0.3.0 below.

## Submission preparation

* Documentation and vignettes are published on the package website
  (<https://alrobles.github.io/xpam/>) rather than shipped in the CRAN
  tarball; `vignettes/` and `inst/doc/` are excluded from the build
  (`.Rbuildignore`), and the `VignetteBuilder` field was dropped so the
  source package carries no prebuilt vignette index.  `README.md` and
  `NEWS.md` are now shipped in the tarball.
* Examples that touch remote data (`get_rodentia_mexico()`,
  `get_world_pca()`, and the downloaded example in `get_geo_pam()`) are
  guarded: the functions document and honour a `"Can't access data."`
  character return when the data server is unreachable, and the examples
  test for it, so `R CMD check --run-donttest` passes offline.  The
  `get_geo_pam()` example additionally runs a fully local `SpatVector`
  example that needs no download at all.
* GDAL remains optional: `configure` detects it and the native raster
  path links only when it is present; without GDAL the package installs,
  loads and the terra-based fallback works.  Engine tests are skipped
  when GDAL is absent.
* Network-dependent tests are gated with `skip_on_cran()` +
  `skip_if_offline()`.

## Engine and API (carried from 0.3.0)

* `pam_env_engine()` streaming environmental PAM; `pam_incidence()`
  sparse incidence (`Matrix::dgCMatrix` or `i`/`j` triplets);
  `pam_jaccard()` species- or site-wise Jaccard with an optional exact
  thresholded edge list; `pam_env()` dense small-input environmental
  path; `has_gdal()` build-time capability report.
* Occurrence sources: in-memory pairs, the `XPOC` binary pair file
  (`xpam_write_occ_bin()`), and GDAL/OGR vector sources -- single paths,
  vectors of paths, directories (k-way merged and deduplicated) or
  in-memory `terra::SpatVector`s.
* OpenMP regions honour `OMP_THREAD_LIMIT`.

# xpam 0.3.0

The engine release: xpam is now a standalone high-performance PAM
package — streamed incidence construction, a sparse PAM surface, and a
native Jaccard backend, on top of the original `get_geo_pam()` API.

## Streaming engine core

* New `pam_env_engine()`: builds the environmental PAM straight from
  raster files. A `PamEngine` C++ object (opaque XPtr) owns GDAL I/O and
  all computation; R passes only paths, the (species, cell) table and
  scalars. Full-width strips are read into a pixel-major buffer, each
  cell is binned by `GridBinning`, and pairs -- pre-sorted by cell --
  merge-join the strips in a single pass into `PamBitset`. Memory is
  O(bitset + strip); `max_bitset_mb` triggers the sort fallback as in
  `pam_env()`. Requires GDAL at build time (`has_gdal()` reports it);
  `configure` degrades gracefully without it. Raster cells are terra
  cellnumbers (1-based).
* `PamBitset` rows are word-aligned (<= 63 pad bits per species), so
  pairwise kernels do `popcount(A & B)/popcount(A | B)` per row pair via
  the `row()` accessor.
* Parallel regions honour `OMP_THREAD_LIMIT` (via `pam::safe_threads()`);
  CRAN check machines set it to 2 and `omp_get_max_threads()` alone does
  not respect it.

## Occurrence sources

* `occ`: data frame of (cell, species) — SEXP-anchored, sorted once
  (`MemOccSource`).
* `occ_file`: cell-sorted "XPOC" binary pair file streamed in ~48 MiB
  chunks (`FileOccSource`); `xpam_write_occ_bin()` writes it. At 4e6
  cells x 500 species the streamed path cuts RSS from ~3.1 GB to ~0.6 GB.
* `occ_polygons`: GDAL/OGR vector sources rasterised in bounded tiles by
  `PolygonOccSource` (`GDALRasterizeGeometries` into a per-feature MEM
  window) — no dense raster, no materialised pair table. Accepts:
  - a path (shapefile, GeoPackage, GeoJSON, ... or an OGR connection
    string such as `duckdb:///db.duckdb` when the GDAL build has the
    driver),
  - a character vector of paths, or a directory scanned for vector
    files — the per-dataset cell-sorted streams are k-way merged and
    deduplicated (`MergeOccSource`); `species` is then required,
  - an in-memory `terra::SpatVector` (bridged through a temporary
    GeoPackage).
  `poly_layer` takes one name per dataset or a single name applied to
  all; `all_touched` selects the burn rule; `occ_dump` tees the stream
  the engine consumes — merged and deduplicated — into a reusable XPOC
  file (`DumpOccSource`).
* Polygon SRS is reprojected to the env raster CRS when they differ;
  features with missing/unmatched fields or unusable geometries are
  counted and reported via `n_skipped` (+ warning).

## Sparse incidence surface

* New `pam_incidence()`: the PAM itself as an R object. Streams any
  occurrence source through the engine and returns the `units x species`
  incidence as a `Matrix::dgCMatrix` (`format = "dgCMatrix"`, Matrix in
  Suggests) or a raw `data.frame(i, j)` triplet table. `space = "geo"`
  is the streaming geographic PAM (raster cells as units, env files as
  grid template only — the scalable replacement for `get_geo_pam()`);
  `space = "env"` gives the environmental incidence over the same bins
  as `pam_env_engine()`. Explicit scale gates: incidence export needs
  the bitmap or sort accumulator (the geo co-occurrence path retains no
  pairs), and `nnz * 16` B must fit `max_out_mb` to marshal to R —
  beyond that the on-disk XPOC from `occ_dump` is the sparse PAM.

## Native Jaccard backend

* `pam_jaccard()`: pairwise Jaccard similarity between species straight
  from any occurrence source — no dense PAM, no RcppParallel.
  `space = "geo"` uses raster cells as incidence units (`env_files`
  supplies only the grid template); `space = "env"` uses env bins.
  Accumulator: the shared bitset when `n_species * n_units/8` fits
  `max_bitset_mb`; in geo space beyond the budget an O(n_species^2)
  co-occurrence matrix is accumulated from per-cell species runs of the
  cell-sorted pair stream (`kernel = "cooc"`, counts returned as
  `cooc`). Semantics port `fastJaccard::jaccard_fast`: J = |A&B| /
  |A|B|, independence expectation `pa*pb/(pa+pb-pa*pb)` on empty unions,
  `centered = TRUE` subtracts the expectation.
* `pam_jaccard(threshold = tau)` returns the similarity graph instead of
  the matrix: an exact prefix-filtered join (AllPairs-style — rare-first
  prefix indexing, size and count bounds, exact merge-join verification)
  emits the `edges` data.frame `(i, j, J >= tau)` with output-sensitive
  work. Same results as filtering the dense matrix; works for both `by`
  modes and all accumulators. On a 2,000-species x 20k-cell fixture it
  is ~80-280x faster than fastJaccard's dense RcppParallel path and
  needs no 305 MB input matrix (`benchmarks/docs/bench_jaccard.md`).
* `pam_jaccard(by = "site")` transposes the question: an
  `n_units x n_units` matrix of similarities between cells or env bins
  over the species they share. Requires the transposed bitset within
  `max_bitset_mb` and the dense `n_units^2` return within `max_out_mb`
  (default 4096 MB, ~22.6k units).

## Package changes

* Package renamed `pam` -> `xpam`: the name `pam` is already taken on
  CRAN. Headers live in `src/*.h` (not `.hpp`) because `R CMD check`
  flags unlikely file names in `src/` as a WARNING.
* Internal refactor: the kernel's heavy state moved to pure-C++ headers
  (`pam_bitset.h`, `pam_binning.h`, `pam_accumulate.h`, `pam_threads.h`)
  with no Rcpp dependency; `SystemRequirements` raised to C++17.
* New `c_bin_assign_d()`: D-dimensional grid binning over a matrix.
* `pam_env()` gains a `max_bitset_mb` memory budget and reports the
  kernel used; accumulation fills a single shared bitset with lock-free
  atomic ORs; out-of-range (species, bin) codes are counted and warned.
* `xpam_sim_env()` / `xpam_sim_ranges()` / `xpam_truth()` /
  `xpam_niche_bound()`: internal deterministic generators and a
  transparent R-only oracle used to validate the engine
  (`test-pam_sim.R`).
* `configure` hardened: `LIBS` is cleared before probing and the GDAL
  link test uses only `gdal-config --libs`, so ambient environment
  libraries (e.g. a user `~/.R/Makevars` exporting `-lRcppParallel` into
  `R CMD INSTALL`'s environment) can no longer silently disable GDAL.
  `AC_PROG_CXX` is also resolved before the link test instead of lazily
  inside it.
* Fixed `url_exists()` returning `TRUE` for unresolvable hosts.
* Benchmarks: `benchmarks/docs/bench_engine.md` (incl. the 1/120-degree global
  sweep: 22.5e9 pairs streamed, ~15-24 GB RSS) and
  `benchmarks/docs/bench_jaccard.md`.

# xpam 0.2.0

* DESCRIPTION overhauled: `Authors@R` format, CRAN-quality title and
  description, `BugReports` URL, `Suggests` for testthat/knitr/rmarkdown.
* Replaced `httr`/`stringr`/`magrittr` dependencies with base R equivalents;
  added `checkmate` for input validation.
* Switched from magrittr pipe (`%>%`) to base R pipe (`|>`).
  Minimum R version raised to 4.1.0.
* Added `checkmate::assert_*` guards in `get_geo_pam()` and `url_exists()`.
* Roxygen documentation rewritten for all exported and internal functions.
* Added testthat test suite (`tests/testthat/`).
* Added introductory vignette (`vignettes/pam-intro.Rmd`).
* Added `NEWS.md`.
* Added R CMD check CI workflow (ubuntu, macOS, windows).
* Fixed `get_rodentia_mexico()` docstring (was incorrectly labelled
  "afrosoricida").

# xpam 0.1.0

* Initial release with `get_geo_pam()`, `get_rodentia_mexico()`, and
  `get_world_pca()`.

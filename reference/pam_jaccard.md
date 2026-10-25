# Pairwise Jaccard similarity between species (native engine)

Streams an occurrence source through the tiled C++ engine and computes
the symmetric `n_species x n_species` Jaccard matrix without
materialising a dense presence-absence matrix: the incidence is held in
a shared bitset (`O(n_species * n_units / 8)` bytes) or – in geographic
space when the bitset would exceed `max_bitset_mb` – in an
`O(n_species^2)` co-occurrence accumulator fed by the cell-sorted pair
stream.

## Usage

``` r
pam_jaccard(
  env_files,
  occ = NULL,
  occ_file = NULL,
  occ_polygons = NULL,
  poly_field = NULL,
  poly_layer = NULL,
  species = NULL,
  all_touched = TRUE,
  occ_dump = NULL,
  space = c("geo", "env"),
  by = c("species", "site"),
  centered = FALSE,
  mins = NULL,
  maxs = NULL,
  nbs = c(30L, 30L),
  poly_chunk_rows = 256L,
  tile_rows = 512,
  threads = 1,
  max_bitset_mb = 1024,
  max_out_mb = 4096,
  threshold = NULL
)
```

## Arguments

- env_files:

  Raster paths (or a file-backed
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)).
  For `space = "geo"` only the grid geometry is used; `occ` cells are
  terra cellnumbers on this grid. For `space = "env"` one multi-band
  file or D single-band files supply the environmental values.

- occ, occ_file, occ_polygons, poly_field, poly_layer, species,
  all_touched, poly_chunk_rows:

  Identical to
  [`pam_env_engine`](https://alrobles.github.io/xpam/reference/pam_env_engine.md).

- occ_dump:

  Optional XPOC dump of the emitted pairs (polygon source only),
  reusable via `occ_file`.

- space:

  `"geo"` (cells; default) or `"env"` (env bins).

- by:

  `"species"` (default; J between species over the units they occupy) or
  `"site"` (J between units over the species they share — the transposed
  view; feasibility-limited, see above).

- centered:

  Subtract the independence expectation (fastJaccard's `center`).
  Default `FALSE`.

- mins, maxs, nbs:

  Binning range and bins per dimension; `space = "env"` only (same
  defaults as
  [`pam_env_engine`](https://alrobles.github.io/xpam/reference/pam_env_engine.md)).

- tile_rows, threads, max_bitset_mb:

  Same as
  [`pam_env_engine`](https://alrobles.github.io/xpam/reference/pam_env_engine.md).

- max_out_mb:

  Feasibility cap (MB) for the dense `n^2` similarity matrix handed back
  to R — independent of `max_bitset_mb` (which bounds the incidence
  state). Default 4096 (~22.6k entities).

- threshold:

  When set to a value in (0, 1\], the engine returns only the pairs with
  similarity \\\ge\\`threshold` as an edge list instead of the dense
  matrix. The join is an exact prefix-filtered similarity join
  (AllPairs-style): entities index only the rarest \\\|A\| - \lceil t_x
  \rceil + 1\\ of their units and candidates are verified exactly, so
  results are identical to filtering the dense matrix — just
  output-sensitive. `NULL` (default) keeps the dense matrix.

## Value

A list with `jaccard` (symmetric matrix — `n_species x n_species` for
`by = "species"`, dimnamed by species; `n_units x n_units` for
`by = "site"`, dimnamed by the unit id: terra cellnumber or env-bin
index, both 1-based), `by`, `species`, `kernel`
(`"bitmap"`/`"sort"`/`"cooc"`), `space`, `nnz`, `footprint` (units per
species), `cooc` (the co-occurrence count matrix, when
`kernel = "cooc"`), `richness` (species per unit, absent under
`"cooc"`), `availability` (cells per env bin, `space = "env"` only),
`n_skipped` and `n_out_of_range`. With `threshold` set, `jaccard` is
replaced by `edges`, a `data.frame(i, j, jaccard)` of the
\\\ge\\`threshold` pairs — `i`/`j` are 1-based entity indices (into
`species` for `by = "species"`, unit ids for `by = "site"`).

## Details

Two incidence spaces are supported. `space = "geo"` uses raster cells as
units (species sets = occupied cells); `env_files` then supplies only
the grid template – its values are never read. With `space = "env"`
units are the D-dimensional environmental bins over `env_files`, exactly
as in
[`pam_env_engine`](https://alrobles.github.io/xpam/reference/pam_env_engine.md).

Values follow fastJaccard semantics: \\J = \|A \cap B\| / \|A \cup
B\|\\; when a pair's union is empty the estimate falls back to the
independence expectation \\E = p_a p_b / (p_a + p_b - p_a p_b)\\ with
\\p = \|A\|/U\\, and `centered = TRUE` subtracts \\E\\. Species
occupying no unit yield `NaN` rows.

`by = "site"` transposes the question: the matrix is then
`n_units x n_units` and `J(i,j)` measures the share of species common to
units \\i\\ and \\j\\. This path always needs the transposed bitset
(units x species) plus the dense `n_units^2` output to fit
`max_bitset_mb` — there is no co-occurrence fallback (it would require a
species-sorted pair stream) — so it is meant for modest grids/binnings;
an informative error is thrown otherwise.

## See also

[`pam_env_engine`](https://alrobles.github.io/xpam/reference/pam_env_engine.md)
(the PAM-construction backend),
[`has_gdal`](https://alrobles.github.io/xpam/reference/has_gdal.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (has_gdal()) {
  # Geographic range-overlap similarity from polygons:
  pam_jaccard("template.tif", occ_polygons = "ranges.gpkg",
              poly_field = "sciname")
}
} # }
```

# Sparse presence-absence incidence from the native engine

Streams an occurrence source through the tiled C++ engine and returns
the incidence matrix itself — the PAM — instead of marginals: a
`n_units x n_species` sparse object built from `(unit, species)` pairs
accumulated in the shared bitset or the sort fallback. This is the
scalable replacement for the dense
[`get_geo_pam`](https://alrobles.github.io/xpam/reference/get_geo_pam.md)
path, which is meant for small inputs.

## Usage

``` r
pam_incidence(
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
  format = c("dgCMatrix", "ij"),
  mins = NULL,
  maxs = NULL,
  nbs = c(30L, 30L),
  poly_chunk_rows = 256L,
  tile_rows = 512,
  threads = 1,
  max_bitset_mb = 1024,
  max_out_mb = 4096
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

- format:

  `"dgCMatrix"` (default; a `units x species` `Matrix::dgCMatrix` of 1s,
  the canonical PAM object) or `"ij"` (a `data.frame(i, j)` of raw
  1-based indices — `i` = unit id: terra cellnumber for `space = "geo"`
  or env-bin index for `space = "env"`; `j` = index into `species`).
  `"dgCMatrix"` requires the Matrix package (in `Suggests`).

- mins, maxs, nbs:

  Binning range and bins per dimension; `space = "env"` only (same
  defaults as
  [`pam_env_engine`](https://alrobles.github.io/xpam/reference/pam_env_engine.md)).

- tile_rows, threads, max_bitset_mb:

  Same as
  [`pam_env_engine`](https://alrobles.github.io/xpam/reference/pam_env_engine.md).

- max_out_mb:

  Feasibility cap (MB) for the returned incidence (nnz \* 16 B as
  triplets). Default 4096 (~268 M pairs). Larger streams belong on disk
  via `occ_dump`.

## Value

A list with `pam` (the sparse incidence object in the requested
`format`), `species`, `kernel` (`"bitmap"`/`"sort"`), `space`, `nnz`,
`richness` (species per unit), `footprint` (units per species), `ncols`,
`nrows`, `n_skipped` and `n_out_of_range`.

## Details

Two incidence spaces are supported. `space = "geo"` uses raster cells as
units (the streaming geographic PAM: each row is a raster cell, i.e. a
terra cellnumber); `env_files` then supplies only the grid template —
its values are never read. With `space = "env"` units are the
D-dimensional environmental bins over `env_files`, exactly as in
[`pam_env_engine`](https://alrobles.github.io/xpam/reference/pam_env_engine.md).

Scale behaviour is explicit rather than silent: the triplets must fit
`max_out_mb` (nnz \* 16 B) to be marshalled to R, the geo co-occurrence
accumulator cannot export incidence (it never retains the pairs), and at
sweep scale (22.5e9 pairs) the honest output is the on-disk XPOC written
by `occ_dump`, not an R object — an informative error is thrown
otherwise.

## See also

[`pam_env_engine`](https://alrobles.github.io/xpam/reference/pam_env_engine.md)
(marginals without the matrix),
[`pam_jaccard`](https://alrobles.github.io/xpam/reference/pam_jaccard.md)
(similarity from the same incidence),
[`get_geo_pam`](https://alrobles.github.io/xpam/reference/get_geo_pam.md)
(dense small-input path),
[`has_gdal`](https://alrobles.github.io/xpam/reference/has_gdal.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (has_gdal()) {
  # Streaming geographic PAM straight from polygons:
  pam_incidence("template.tif", occ_polygons = "ranges.gpkg",
                poly_field = "sciname")
}
} # }
```

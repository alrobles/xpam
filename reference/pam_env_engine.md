# Build a sparse environmental PAM directly from raster files (GDAL engine)

Runs the tiled C++ engine: environmental rasters are read strip by strip
with GDAL, each cell is assigned to a climate bin, and the (species,
cell) pairs are scattered into a shared bitset – all inside C++. R only
passes file paths, the pair table and scalars; the raster layers are
never loaded into the R session, and no dense species-by-site or
pair-by-bin object of length P is materialised.

## Usage

``` r
pam_env_engine(
  env_files,
  occ = NULL,
  occ_file = NULL,
  occ_polygons = NULL,
  poly_field = NULL,
  poly_layer = NULL,
  species = NULL,
  all_touched = TRUE,
  occ_dump = NULL,
  poly_chunk_rows = 256L,
  mins = NULL,
  maxs = NULL,
  nbs,
  tile_rows = 512,
  threads = 1,
  max_bitset_mb = 1024
)
```

## Arguments

- env_files:

  Character vector of raster paths: either one multi-band file (D =
  number of bands) or D single-band files. Alternatively a
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  (sources are extracted). All bands must share raster dimensions.

- occ:

  A data frame with columns `cell` and `species`, where `cell` is the
  1-based raster cell index (as in
  [`cells`](https://rspatial.github.io/terra/reference/cells.html) /
  [`cellFromXY`](https://rspatial.github.io/terra/reference/xyCellFrom.html)).
  Mutually exclusive with `occ_file`.

- occ_file:

  Path to an XPOC occurrence file: a cell-sorted binary pair table of
  (cell, species) records, such as the one produced by the `occ_dump`
  argument. When given, the pair table is streamed from disk in bounded
  chunks and never materialised in R or C++.

- occ_polygons:

  Vector dataset(s) readable by GDAL — a path to a shapefile,
  GeoPackage, GeoJSON or any other OGR driver (an OGR connection string
  such as `duckdb:///db.duckdb` works when the GDAL build has the
  driver); a character vector of such paths; a directory, which is
  scanned for vector files; or an in-memory
  [`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html)
  (bridged through a temporary GeoPackage, e.g. the result of a
  DuckDB/sf query). With more than one dataset the cell-sorted streams
  are k-way merged and deduplicated, and `species` is then required
  (levels cannot be derived consistently across files). Geometries are
  rasterised onto the `env_files` grid chunk by chunk and streamed as
  (species, cell) pairs; no dense intermediate is materialised. Mutually
  exclusive with `occ` and `occ_file`.

- poly_field:

  Attribute of `occ_polygons` identifying the species — the same field
  name is used in every dataset. Values are read as strings and matched
  against `species`; required when `occ_polygons` is given.

- poly_layer:

  Layer name(s) inside `occ_polygons`: either a single name applied to
  every dataset, or one name per dataset; the default `NULL` uses the
  first layer of each.

- species:

  Character vector giving the species order (code = index). With a
  single `occ_polygons` dataset, `NULL` derives the levels from the
  field's first occurrence in the file; with several datasets `species`
  is required. Values not in `species` are skipped (and counted in
  `n_skipped`).

- all_touched:

  Logical: burn every pixel a polygon touches (`TRUE`, the GDAL
  `ALL_TOUCHED` rule) or only pixels whose centre lies inside (`FALSE`).
  Default `TRUE`.

- occ_dump:

  Optional path where the emitted (species, cell) pairs are written as
  an XPOC file while streaming – the reproducibility artefact, reusable
  later via `occ_file`.

- poly_chunk_rows:

  Rows per rasterise/emit chunk inside the polygon source. Default
  `256`.

- mins, maxs:

  Numeric vectors of length D: binning range per environmental
  dimension. If `NULL`, the 0.5%/99.5% quantiles of each raster band are
  used (requires terra).

- nbs:

  Integer vector of length D: number of bins per dimension. For `nb` =
  scalar shortcut, recycle it to length D.

- tile_rows:

  Integer: rows per strip. Default `512`.

- threads:

  Integer: OpenMP threads for binning/scatter (honours
  `OMP_THREAD_LIMIT`). Default `1`.

- max_bitset_mb:

  Numeric: bitset memory budget in MiB before the sort kernel is used.
  Default `1024`.

## Value

A list with `nnz`, `richness` (distinct species per bin), `footprint`
(bins per species), `availability` (cells per bin), `kernel`,
`state_mb`, `species`, `nbs`, `ncols`, `nrows` and `n_strips`; plus the
counters `n_skipped` (polygon features dropped: missing field, unmatched
species or unusable geometry) and `n_out_of_range` (pairs outside the
binning extent).

## Details

Memory is `O(bitset + strip)`: the bitset is `n_species * n_bins / 8`
bytes and the strip is `tile_rows * ncols * D` doubles. When the bitset
would exceed `max_bitset_mb`, the engine falls back to a sort-based
kernel (memory `O(P)`) and reports which kernel ran.

## See also

[`pam_env`](https://alrobles.github.io/xpam/reference/pam_env.md)
(in-memory path),
[`has_gdal`](https://alrobles.github.io/xpam/reference/has_gdal.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (has_gdal()) {
  res <- pam_env_engine("clim_2band.tif", occ, nbs = c(30, 30))
}
} # }
```

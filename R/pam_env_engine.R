#' Build a sparse environmental PAM directly from raster files (GDAL engine)
#'
#' Runs the tiled C++ engine: environmental rasters are read strip by strip
#' with GDAL, each cell is assigned to a climate bin, and the (species, cell)
#' pairs are scattered into a shared bitset -- all inside C++.  R only passes
#' file paths, the pair table and scalars; the raster layers are never loaded
#' into the R session, and no dense species-by-site or pair-by-bin object of
#' length P is materialised.
#'
#' Memory is \code{O(bitset + strip)}: the bitset is
#' \code{n_species * n_bins / 8} bytes and the strip is
#' \code{tile_rows * ncols * D} doubles.  When the bitset would exceed
#' \code{max_bitset_mb}, the engine falls back to a sort-based kernel
#' (memory \code{O(P)}) and reports which kernel ran.
#'
#' @param env_files Character vector of raster paths: either one multi-band
#'   file (D = number of bands) or D single-band files.  Alternatively a
#'   \code{\link[terra]{SpatRaster}} (sources are extracted).  All bands must
#'   share raster dimensions.
#' @param occ A data frame with columns \code{cell} and \code{species}, where
#'   \code{cell} is the 1-based raster cell index (as in
#'   \code{\link[terra]{cells}} / \code{\link[terra]{cellFromXY}}).  Mutually
#'   exclusive with \code{occ_file}.
#' @param occ_file Path to an XPOC occurrence file: a cell-sorted binary
#'   pair table of (cell, species) records, such as the one produced by the
#'   \code{occ_dump} argument.  When given, the pair table is streamed from
#'   disk in bounded chunks and never materialised in R or C++.
#' @param occ_polygons Vector dataset(s) readable by GDAL — a path to a
#'   shapefile, GeoPackage, GeoJSON or any other OGR driver (an OGR
#'   connection string such as \code{duckdb:///db.duckdb} works when the
#'   GDAL build has the driver); a character vector of such paths; a
#'   directory, which is scanned for vector files; or an in-memory
#'   \code{terra::SpatVector} (bridged through a temporary GeoPackage,
#'   e.g. the result of a DuckDB/\pkg{sf} query).  With more than one
#'   dataset the cell-sorted streams are k-way merged and deduplicated,
#'   and \code{species} is then required (levels cannot be derived
#'   consistently across files).  Geometries are rasterised onto the
#'   \code{env_files} grid chunk by chunk and streamed as (species, cell)
#'   pairs; no dense intermediate is materialised.  Mutually exclusive
#'   with \code{occ} and \code{occ_file}.
#' @param poly_field Attribute of \code{occ_polygons} identifying the
#'   species — the same field name is used in every dataset.  Values are
#'   read as strings and matched against \code{species}; required when
#'   \code{occ_polygons} is given.
#' @param poly_layer Layer name(s) inside \code{occ_polygons}: either a
#'   single name applied to every dataset, or one name per dataset; the
#'   default \code{NULL} uses the first layer of each.
#' @param species Character vector giving the species order (code = index).
#'   With a single \code{occ_polygons} dataset, \code{NULL} derives the
#'   levels from the field's first occurrence in the file; with several
#'   datasets \code{species} is required.  Values not in \code{species}
#'   are skipped (and counted in \code{n_skipped}).
#' @param all_touched Logical: burn every pixel a polygon touches
#'   (\code{TRUE}, the GDAL \code{ALL_TOUCHED} rule) or only pixels whose
#'   centre lies inside (\code{FALSE}).  Default \code{TRUE}.
#' @param occ_dump Optional path where the emitted (species, cell) pairs
#'   are written as an XPOC file while streaming -- the reproducibility
#'   artefact, reusable later via \code{occ_file}.
#' @param poly_chunk_rows Rows per rasterise/emit chunk inside the polygon
#'   source.  Default \code{256}.
#' @param mins,maxs Numeric vectors of length D: binning range per
#'   environmental dimension.  If \code{NULL}, the 0.5\%/99.5\% quantiles of
#'   each raster band are used (requires \pkg{terra}).
#' @param nbs Integer vector of length D: number of bins per dimension.  For
#'   \code{nb} = scalar shortcut, recycle it to length D.
#' @param tile_rows Integer: rows per strip.  Default \code{512}.
#' @param threads Integer: OpenMP threads for binning/scatter (honours
#'   \code{OMP_THREAD_LIMIT}).  Default \code{1}.
#' @param max_bitset_mb Numeric: bitset memory budget in MiB before the sort
#'   kernel is used.  Default \code{1024}.
#'
#' @return A list with \code{nnz}, \code{richness} (distinct species per bin),
#'   \code{footprint} (bins per species), \code{availability} (cells per
#'   bin), \code{kernel}, \code{state_mb}, \code{species}, \code{nbs},
#'   \code{ncols}, \code{nrows} and \code{n_strips}; plus the counters
#'   \code{n_skipped} (polygon features dropped: missing field, unmatched
#'   species or unusable geometry) and \code{n_out_of_range} (pairs outside
#'   the binning extent).
#'
#' @seealso \code{\link{pam_env}} (in-memory path), \code{\link{has_gdal}}
#' @export
#'
#' @examples
#' \dontrun{
#' if (has_gdal()) {
#'   res <- pam_env_engine("clim_2band.tif", occ, nbs = c(30, 30))
#' }
#' }
pam_env_engine <- function(env_files, occ = NULL, occ_file = NULL,
                           occ_polygons = NULL, poly_field = NULL,
                           poly_layer = NULL, species = NULL,
                           all_touched = TRUE, occ_dump = NULL,
                           poly_chunk_rows = 256L,
                           mins = NULL, maxs = NULL, nbs,
                           tile_rows = 512, threads = 1, max_bitset_mb = 1024) {
  if (!has_gdal())
    stop("xpam was built without GDAL support; use pam_env() instead",
         call. = FALSE)
  env_files <- .pam_env_paths(env_files)
  checkmate::assert_integerish(nbs, min.len = 1, lower = 1)

  # Per-dim quantile defaults need band summaries; terra reads lazily.
  if (is.null(mins) || is.null(maxs)) {
    r <- terra::rast(env_files)
    qs <- terra::global(r, fun = function(v)
                        stats::quantile(v, probs = c(0.005, 0.995),
                                        na.rm = TRUE))
    mins <- qs[[1]]; maxs <- qs[[2]]
  }
  checkmate::assert_numeric(mins, len = length(nbs), finite = TRUE)
  checkmate::assert_numeric(maxs, len = length(nbs), finite = TRUE)
  checkmate::assert_number(tile_rows, lower = 1)
  checkmate::assert_number(threads, lower = 1)
  checkmate::assert_number(max_bitset_mb, lower = 0)

  eng <- xpam_engine_create()
  xpam_engine_open_env(eng, env_files)
  species_lv <- .pam_engine_source(eng, occ, occ_file, occ_polygons,
                                   poly_field, poly_layer, species,
                                   all_touched, occ_dump, poly_chunk_rows)
  xpam_engine_set_bins(eng, as.numeric(mins), as.numeric(maxs),
                       as.integer(nbs))
  xpam_engine_set_tile_rows(eng, as.integer(tile_rows))
  xpam_engine_set_threads(eng, as.integer(threads))
  xpam_engine_set_budget_mb(eng, as.numeric(max_bitset_mb))
  res <- xpam_engine_compute(eng)

  if (res$n_out_of_range > 0)
    warning(sprintf("%d (species, cell) pair(s) were out of range and were skipped",
                    res$n_out_of_range))
  if (res$n_skipped > 0)
    warning(sprintf("%d polygon feature(s) were skipped (missing field, unmatched species or unusable geometry)",
                    res$n_skipped))
  if (is.null(species_lv) && length(res$species_levels))
    species_lv <- res$species_levels
  res$species_levels <- NULL
  res$species <- species_lv
  res$nbs <- as.integer(nbs)
  res
}

#' Check whether the package was built with GDAL support
#'
#' @return Logical scalar: \code{TRUE} when the native GDAL engine
#'   (\code{\link{pam_env_engine}}) is available.
#' @export
#' @examples
#' has_gdal()
has_gdal <- function() {
  isTRUE(xpam_have_gdal())
}

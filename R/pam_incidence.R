#' Sparse presence-absence incidence from the native engine
#'
#' Streams an occurrence source through the tiled C++ engine and returns
#' the incidence matrix itself — the PAM — instead of marginals: a
#' \code{n_units x n_species} sparse object built from
#' \code{(unit, species)} pairs accumulated in the shared bitset or the
#' sort fallback.  This is the scalable replacement for the dense
#' \code{\link{get_geo_pam}} path, which is meant for small inputs.
#'
#' Two incidence spaces are supported.  \code{space = "geo"} uses raster
#' cells as units (the streaming geographic PAM: each row is a raster
#' cell, i.e. a terra cellnumber); \code{env_files} then supplies only
#' the grid template — its values are never read.  With
#' \code{space = "env"} units are the D-dimensional environmental bins
#' over \code{env_files}, exactly as in \code{\link{pam_env_engine}}.
#'
#' Scale behaviour is explicit rather than silent: the triplets must fit
#' \code{max_out_mb} (nnz * 16 B) to be marshalled to R, the
#' geo co-occurrence accumulator cannot export incidence (it never
#' retains the pairs), and at sweep scale (22.5e9 pairs) the honest
#' output is the on-disk XPOC written by \code{occ_dump}, not an R
#' object — an informative error is thrown otherwise.
#'
#' @param env_files Raster paths (or a file-backed
#'   \code{\link[terra]{SpatRaster}}).  For \code{space = "geo"} only the
#'   grid geometry is used; \code{occ} cells are terra cellnumbers on
#'   this grid.  For \code{space = "env"} one multi-band file or D
#'   single-band files supply the environmental values.
#' @param occ,occ_file,occ_polygons,poly_field,poly_layer,species,all_touched,poly_chunk_rows
#'   Identical to \code{\link{pam_env_engine}}.
#' @param occ_dump Optional XPOC dump of the emitted pairs (polygon
#'   source only), reusable via \code{occ_file}.
#' @param space \code{"geo"} (cells; default) or \code{"env"} (env bins).
#' @param format \code{"dgCMatrix"} (default; a \code{units x species}
#'   \code{Matrix::dgCMatrix} of 1s, the canonical PAM object) or
#'   \code{"ij"} (a \code{data.frame(i, j)} of raw 1-based indices —
#'   \code{i} = unit id: terra cellnumber for \code{space = "geo"} or
#'   env-bin index for \code{space = "env"}; \code{j} = index into
#'   \code{species}).  \code{"dgCMatrix"} requires the \pkg{Matrix}
#'   package (in \code{Suggests}).
#' @param mins,maxs,nbs Binning range and bins per dimension;
#'   \code{space = "env"} only (same defaults as
#'   \code{\link{pam_env_engine}}).
#' @param tile_rows,threads,max_bitset_mb Same as
#'   \code{\link{pam_env_engine}}.
#' @param max_out_mb Feasibility cap (MB) for the returned incidence
#'   (nnz * 16 B as triplets).  Default 4096 (~268 M pairs).  Larger
#'   streams belong on disk via \code{occ_dump}.
#'
#' @return A list with \code{pam} (the sparse incidence object in the
#'   requested \code{format}), \code{species}, \code{kernel}
#'   (\code{"bitmap"}/\code{"sort"}), \code{space}, \code{nnz},
#'   \code{richness} (species per unit), \code{footprint} (units per
#'   species), \code{ncols}, \code{nrows}, \code{n_skipped} and
#'   \code{n_out_of_range}.
#'
#' @seealso \code{\link{pam_env_engine}} (marginals without the matrix),
#'   \code{\link{pam_jaccard}} (similarity from the same incidence),
#'   \code{\link{get_geo_pam}} (dense small-input path), \code{\link{has_gdal}}
#' @export
#'
#' @examples
#' \dontrun{
#' if (has_gdal()) {
#'   # Streaming geographic PAM straight from polygons:
#'   pam_incidence("template.tif", occ_polygons = "ranges.gpkg",
#'                 poly_field = "sciname")
#' }
#' }
pam_incidence <- function(env_files, occ = NULL, occ_file = NULL,
                          occ_polygons = NULL, poly_field = NULL,
                          poly_layer = NULL, species = NULL,
                          all_touched = TRUE, occ_dump = NULL,
                          space = c("geo", "env"),
                          format = c("dgCMatrix", "ij"),
                          mins = NULL, maxs = NULL, nbs = c(30L, 30L),
                          poly_chunk_rows = 256L,
                          tile_rows = 512, threads = 1,
                          max_bitset_mb = 1024, max_out_mb = 4096) {
  if (!has_gdal())
    stop("xpam was built without GDAL support", call. = FALSE)
  space <- match.arg(space)
  format <- match.arg(format)
  env_files <- .pam_env_paths(env_files)
  checkmate::assert_integerish(nbs, min.len = 1, lower = 1)

  if (format == "dgCMatrix" && !requireNamespace("Matrix", quietly = TRUE))
    stop("format = 'dgCMatrix' needs the Matrix package; ",
         "install it or use format = 'ij'", call. = FALSE)

  if (space == "env" && (is.null(mins) || is.null(maxs))) {
    r <- terra::rast(env_files)
    qs <- terra::global(r, fun = function(v)
                        stats::quantile(v, probs = c(0.005, 0.995),
                                        na.rm = TRUE))
    mins <- qs[[1]]; maxs <- qs[[2]]
  }
  if (space == "env") {
    checkmate::assert_numeric(mins, len = length(nbs), finite = TRUE)
    checkmate::assert_numeric(maxs, len = length(nbs), finite = TRUE)
  }
  checkmate::assert_number(tile_rows, lower = 1)
  checkmate::assert_number(threads, lower = 1)
  checkmate::assert_number(max_bitset_mb, lower = 0)
  checkmate::assert_number(max_out_mb, lower = 0)

  eng <- xpam_engine_create()
  xpam_engine_open_env(eng, env_files)
  species_lv <- .pam_engine_source(eng, occ, occ_file, occ_polygons,
                                   poly_field, poly_layer, species,
                                   all_touched, occ_dump, poly_chunk_rows)
  if (space == "env")
    xpam_engine_set_bins(eng, as.numeric(mins), as.numeric(maxs),
                         as.integer(nbs))
  xpam_engine_set_space(eng, if (space == "geo") 1L else 0L)
  xpam_engine_set_export_incidence(eng, TRUE)
  xpam_engine_set_tile_rows(eng, as.integer(tile_rows))
  xpam_engine_set_threads(eng, as.integer(threads))
  xpam_engine_set_budget_mb(eng, as.numeric(max_bitset_mb))
  xpam_engine_set_out_budget_mb(eng, as.numeric(max_out_mb))
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

  n_units <- if (space == "geo")
    as.numeric(res$ncols) * as.numeric(res$nrows)
  else
    prod(nbs)

  # [[ ]] guards against partial matching of optional elements
  ii <- res[["inc_i"]]; jj <- res[["inc_j"]]
  if (format == "dgCMatrix") {
    res$pam <- Matrix::sparseMatrix(
      i = ii, j = jj, x = rep(1, length(ii)),
      dims = c(n_units, length(species_lv)),
      dimnames = list(NULL, species_lv))
  } else {
    res$pam <- data.frame(i = ii, j = jj)
  }
  res$inc_i <- res$inc_j <- NULL
  res
}

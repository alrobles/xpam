#' Configure an occurrence source on a PamEngine instance
#'
#' Shared validation + engine calls for `pam_env_engine()` and
#' `pam_jaccard()`: exactly one of `occ`/`occ_file`/`occ_polygons`.
#'
#' @return The species-level vector in code order (`NULL` when it will be
#'   derived inside C++), invisibly the engine is configured.
#' @noRd
.pam_engine_source <- function(eng, occ = NULL, occ_file = NULL,
                               occ_polygons = NULL, poly_field = NULL,
                               poly_layer = NULL, species = NULL,
                               all_touched = TRUE, occ_dump = NULL,
                               poly_chunk_rows = 256L) {
  n_src <- sum(!is.null(occ), !is.null(occ_file), !is.null(occ_polygons))
  if (n_src != 1L)
    stop("pass exactly one of `occ`, `occ_file` or `occ_polygons`",
         call. = FALSE)
  species_lv <- NULL
  if (!is.null(occ)) {
    checkmate::assert_data_frame(occ)
    stopifnot(all(c("cell", "species") %in% names(occ)))
    species_lv <- unique(occ$species)
    spp <- match(occ$species, species_lv) - 1L
    cell0 <- as.integer(occ$cell) - 1L        # terra cellnumber -> 0-based
    xpam_engine_set_occ(eng, spp, cell0, length(species_lv))
  } else if (!is.null(occ_file)) {
    checkmate::assert_file_exists(occ_file)
    xpam_engine_set_occ_file(eng, occ_file)   # levels live in the header
  } else {
    occ_polygons <- .pam_poly_paths(occ_polygons)
    n_ds <- length(occ_polygons)
    checkmate::assert_string(poly_field)
    checkmate::assert_character(poly_layer, null.ok = TRUE, min.len = 1,
                                max.len = n_ds, any.missing = FALSE)
    if (!is.null(poly_layer) && length(poly_layer) != n_ds)
      poly_layer <- rep_len(poly_layer, n_ds)   # scalar -> recycled
    checkmate::assert_character(species, null.ok = TRUE, min.len = 1,
                                any.missing = FALSE, unique = TRUE)
    if (n_ds > 1 && is.null(species))
      stop("`species` must be given when `occ_polygons` has multiple ",
           "datasets (levels cannot be derived consistently across files)",
           call. = FALSE)
    checkmate::assert_flag(all_touched)
    if (!is.null(occ_dump))
      checkmate::assert_path_for_output(occ_dump, overwrite = TRUE)
    checkmate::assert_count(poly_chunk_rows, positive = TRUE)
    xpam_engine_set_occ_polygons(
      eng, occ_polygons,
      if (is.null(poly_layer)) character(0) else poly_layer,
      poly_field, if (is.null(species)) character(0) else species,
      all_touched, as.integer(poly_chunk_rows),
      if (is.null(occ_dump)) "" else occ_dump)
    species_lv <- species                     # NULL -> levels from C++
  }
  species_lv
}

#' Resolve polygon inputs to OGR dataset paths
#'
#' `occ_polygons` accepts: a vector-dataset path (any OGR driver — .shp,
#' .gpkg, .geojson, ... or a DSN like `duckdb:///db.duckdb` when GDAL has
#' the driver), a character vector of such paths (merged by the engine),
#' a directory scanned for vector files, or an in-memory
#' `terra::SpatVector` (bridged through a temporary .gpkg).
#'
#' @return Character vector of dataset paths.
#' @noRd
.pam_poly_paths <- function(occ_polygons) {
  if (inherits(occ_polygons, "SpatVector")) {
    f <- tempfile(fileext = ".gpkg")
    terra::writeVector(occ_polygons, f, filetype = "GPKG",
                       overwrite = TRUE)
    return(f)
  }
  checkmate::assert_character(occ_polygons, min.len = 1,
                              any.missing = FALSE)
  if (length(occ_polygons) == 1 && dir.exists(occ_polygons)) {
    paths <- list.files(
      occ_polygons,
      pattern = "\\.(shp|gpkg|geojson|json|fgb|sqlite|gml|kml|gpkg2)$",
      full.names = TRUE, ignore.case = TRUE)
    paths <- paths[!dir.exists(paths)]
    if (!length(paths))
      stop("no OGR vector files found in directory: ", occ_polygons,
           call. = FALSE)
    return(sort(paths))
  }
  for (p in occ_polygons) {
    # OGR connection strings (duckdb:///, PG:, /vsizip/...) need no file
    if (grepl("://", p, fixed = TRUE) || grepl("^/vsi", p)) next
    if (grepl("^[A-Za-z]+:", p) && !file.exists(p)) next  # e.g. PG:
    checkmate::assert_file_exists(p)
  }
  occ_polygons
}

#' Resolve raster-file paths (vector of paths or a file-backed SpatRaster)
#' @noRd
.pam_env_paths <- function(env_files) {
  if (inherits(env_files, "SpatRaster")) {
    src <- terra::sources(env_files)
    if (any(!nzchar(src) | !file.exists(src)))
      stop("SpatRaster is not file-backed; write it to disk first",
           call. = FALSE)
    env_files <- src
  }
  checkmate::assert_character(env_files, min.len = 1, any.missing = FALSE)
  if (any(!file.exists(env_files)))
    stop("env file not found: ", env_files[!file.exists(env_files)][1],
         call. = FALSE)
  env_files
}

#' Pairwise Jaccard similarity between species (native engine)
#'
#' Streams an occurrence source through the tiled C++ engine and computes
#' the symmetric \code{n_species x n_species} Jaccard matrix without
#' materialising a dense presence-absence matrix: the incidence is held in
#' a shared bitset (\code{O(n_species * n_units / 8)} bytes) or -- in
#' geographic space when the bitset would exceed \code{max_bitset_mb} -- in
#' an \code{O(n_species^2)} co-occurrence accumulator fed by the
#' cell-sorted pair stream.
#'
#' Two incidence spaces are supported.  \code{space = "geo"} uses raster
#' cells as units (species sets = occupied cells); \code{env_files} then
#' supplies only the grid template -- its values are never read.  With
#' \code{space = "env"} units are the D-dimensional environmental bins over
#' \code{env_files}, exactly as in \code{\link{pam_env_engine}}.
#'
#' Values follow \pkg{fastJaccard} semantics: \eqn{J = |A \cap B| /
#' |A \cup B|}; when a pair's union is empty the estimate falls back to the
#' independence expectation \eqn{E = p_a p_b / (p_a + p_b - p_a p_b)} with
#' \eqn{p = |A|/U}, and \code{centered = TRUE} subtracts \eqn{E}.  Species
#' occupying no unit yield \code{NaN} rows.
#'
#' \code{by = "site"} transposes the question: the matrix is then
#' \code{n_units x n_units} and \code{J(i,j)} measures the share of
#' species common to units \eqn{i} and \eqn{j}.  This path always needs
#' the transposed bitset (units x species) plus the dense
#' \code{n_units^2} output to fit \code{max_bitset_mb} — there is no
#' co-occurrence fallback (it would require a species-sorted pair
#' stream) — so it is meant for modest grids/binnings; an informative
#' error is thrown otherwise.
#'
#' @param env_files Raster paths (or a file-backed
#'   \code{\link[terra]{SpatRaster}}).  For \code{space = "geo"} only the
#'   grid geometry is used; \code{occ} cells are terra cellnumbers on this
#'   grid.  For \code{space = "env"} one multi-band file or D single-band
#'   files supply the environmental values.
#' @param occ,occ_file,occ_polygons,poly_field,poly_layer,species,all_touched,poly_chunk_rows
#'   Identical to \code{\link{pam_env_engine}}.
#' @param occ_dump Optional XPOC dump of the emitted pairs (polygon source
#'   only), reusable via \code{occ_file}.
#' @param space \code{"geo"} (cells; default) or \code{"env"} (env bins).
#' @param by \code{"species"} (default; J between species over the units
#'   they occupy) or \code{"site"} (J between units over the species they
#'   share — the transposed view; feasibility-limited, see above).
#' @param centered Subtract the independence expectation (fastJaccard's
#'   \code{center}).  Default \code{FALSE}.
#' @param mins,maxs,nbs Binning range and bins per dimension;
#'   \code{space = "env"} only (same defaults as
#'   \code{\link{pam_env_engine}}).
#' @param tile_rows,threads,max_bitset_mb Same as
#'   \code{\link{pam_env_engine}}.
#' @param max_out_mb Feasibility cap (MB) for the dense \code{n^2}
#'   similarity matrix handed back to R — independent of
#'   \code{max_bitset_mb} (which bounds the incidence state).  Default
#'   4096 (~22.6k entities).
#' @param threshold When set to a value in (0, 1], the engine returns
#'   only the pairs with similarity \eqn{\ge}\code{threshold} as an edge
#'   list instead of the dense matrix.  The join is an exact
#'   prefix-filtered similarity join (AllPairs-style): entities index
#'   only the rarest \eqn{|A| - \lceil t_x \rceil + 1} of their units and
#'   candidates are verified exactly, so results are identical to
#'   filtering the dense matrix — just output-sensitive.  \code{NULL}
#'   (default) keeps the dense matrix.
#'
#' @return A list with \code{jaccard} (symmetric matrix —
#'   \code{n_species x n_species} for \code{by = "species"}, dimnamed by
#'   species; \code{n_units x n_units} for \code{by = "site"}, dimnamed
#'   by the unit id: terra cellnumber or env-bin index, both 1-based),
#'   \code{by}, \code{species}, \code{kernel}
#'   (\code{"bitmap"}/\code{"sort"}/\code{"cooc"}), \code{space},
#'   \code{nnz}, \code{footprint} (units per species), \code{cooc} (the
#'   co-occurrence count matrix, when \code{kernel = "cooc"}),
#'   \code{richness} (species per unit, absent under \code{"cooc"}),
#'   \code{availability} (cells per env bin, \code{space = "env"} only),
#'   \code{n_skipped} and \code{n_out_of_range}.  With
#'   \code{threshold} set, \code{jaccard} is replaced by \code{edges}, a
#'   \code{data.frame(i, j, jaccard)} of the \eqn{\ge}\code{threshold}
#'   pairs — \code{i}/\code{j} are 1-based entity indices (into
#'   \code{species} for \code{by = "species"}, unit ids for
#'   \code{by = "site"}).
#'
#' @seealso \code{\link{pam_env_engine}} (the PAM-construction backend),
#'   \code{\link{has_gdal}}
#' @export
#'
#' @examples
#' \dontrun{
#' if (has_gdal()) {
#'   # Geographic range-overlap similarity from polygons:
#'   pam_jaccard("template.tif", occ_polygons = "ranges.gpkg",
#'               poly_field = "sciname")
#' }
#' }
pam_jaccard <- function(env_files, occ = NULL, occ_file = NULL,
                        occ_polygons = NULL, poly_field = NULL,
                        poly_layer = NULL, species = NULL,
                        all_touched = TRUE, occ_dump = NULL,
                        space = c("geo", "env"), by = c("species", "site"),
                        centered = FALSE,
                        mins = NULL, maxs = NULL, nbs = c(30L, 30L),
                        poly_chunk_rows = 256L,
                        tile_rows = 512, threads = 1, max_bitset_mb = 1024,
                        max_out_mb = 4096, threshold = NULL) {
  if (!has_gdal())
    stop("xpam was built without GDAL support", call. = FALSE)
  space <- match.arg(space)
  by <- match.arg(by)
  env_files <- .pam_env_paths(env_files)
  checkmate::assert_flag(centered)
  checkmate::assert_integerish(nbs, min.len = 1, lower = 1)

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
  checkmate::assert_number(threshold, lower = 0, upper = 1, null.ok = TRUE)

  eng <- xpam_engine_create()
  xpam_engine_open_env(eng, env_files)
  species_lv <- .pam_engine_source(eng, occ, occ_file, occ_polygons,
                                   poly_field, poly_layer, species,
                                   all_touched, occ_dump, poly_chunk_rows)
  if (space == "env")
    xpam_engine_set_bins(eng, as.numeric(mins), as.numeric(maxs),
                         as.integer(nbs))
  xpam_engine_set_space(eng, if (space == "geo") 1L else 0L)
  xpam_engine_set_jaccard(eng, if (centered) 2L else 1L)
  xpam_engine_set_by(eng, if (by == "site") 1L else 0L)
  xpam_engine_set_tile_rows(eng, as.integer(tile_rows))
  xpam_engine_set_threads(eng, as.integer(threads))
  xpam_engine_set_budget_mb(eng, as.numeric(max_bitset_mb))
  xpam_engine_set_out_budget_mb(eng, as.numeric(max_out_mb))
  if (!is.null(threshold))
    xpam_engine_set_jaccard_threshold(eng, as.numeric(threshold))
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
  # NB: [[ ]] everywhere below -- res$jaccard would partial-match
  # jaccard_mode (and res$species, species_levels) when the element is
  # absent.
  if (!is.null(res[["jaccard"]])) {
    n_ent <- nrow(res$jaccard)
    dimnames(res$jaccard) <- if (by == "species" && !is.null(species_lv))
      list(species_lv, species_lv)
    else
      list(seq_len(n_ent), seq_len(n_ent))  # unit ids (cellnumber/bin)
  }
  if (!is.null(res[["edge_i"]])) {
    res$edges <- data.frame(i = res$edge_i, j = res$edge_j,
                            jaccard = res$edge_v)
    res$edge_i <- res$edge_j <- res$edge_v <- NULL
  }
  if (!is.null(res[["cooc"]]) && !is.null(species_lv))
    dimnames(res$cooc) <- list(species_lv, species_lv)
  res
}

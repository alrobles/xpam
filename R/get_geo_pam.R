#' Build a Presence-Absence Matrix from Species Range Polygons
#'
#' Rasterizes a set of species range polygons onto a regular grid and
#' returns a sites-by-species presence-absence matrix (PAM) as a tibble.
#'
#' @param shp A \code{\link[terra]{SpatVector}} with a \code{sciname}
#'   attribute column identifying each species.
#' @param res Numeric.
#'   Grid cell size in the same units as the coordinate reference system
#'   of \code{shp} (typically decimal degrees). Default \code{0.5}.
#' @param ... Additional arguments passed to \code{\link[terra]{rasterize}}.
#'
#' @return A \code{\link[tibble]{tibble}} with columns \code{x}, \code{y}
#'   (cell centroids) followed by one column per species containing
#'   \code{1} (present) or \code{0} (absent).
#'
#' @export
#'
#' @examples
#' \donttest{
#' # A small, fully local example: two square ranges, no download needed
#' square <- function(x0, y0, size = 2) {
#'   m <- matrix(c(x0, y0, x0 + size, y0,
#'                 x0 + size, y0 + size, x0, y0 + size),
#'               ncol = 2, byrow = TRUE)
#'   terra::vect(m, type = "polygons", crs = "EPSG:4326")
#' }
#' ranges <- rbind(square(0, 0), square(2, 2))
#' ranges$sciname <- c("species_a", "species_b")
#' pam <- get_geo_pam(ranges, res = 1)
#' head(pam)
#'
#' # The packaged example data are downloaded on demand, so guard for
#' # machines without network access (the function returns a character
#' # string instead of a SpatVector when the server is unreachable)
#' shp <- get_rodentia_mexico()
#' if (!is.character(shp)) {
#'   pam_mex <- get_geo_pam(shp, res = 1)
#'   head(pam_mex)
#' }
#' }
get_geo_pam <- function(shp, res = 0.5, ...) {
  checkmate::assert_class(shp, "SpatVector")
  checkmate::assert_number(res, lower = 0, finite = TRUE)

  r <- terra::rast(shp, res = res)
  z <- terra::rasterize(shp, r, "sciname")
  z_df <- terra::as.data.frame(z, xy = TRUE)

  z_df |>
    dplyr::mutate(value = 1) |>
    tidyr::spread(key = 3, value = 4, fill = 0) |>
    tibble::as_tibble()
}

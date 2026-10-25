#' Retrieve a WorldClim PCA Environmental Raster
#'
#' Downloads a \code{\link[terra]{SpatRaster}} containing the first two
#' principal components derived from WorldClim 2.0 bioclimatic variables
#' at global extent.
#'
#' @return A \code{\link[terra]{SpatRaster}} object, or
#'   \code{"Can't access data."} if the remote server is unreachable.
#'
#' @export
#'
#' @source WorldClim 2.0 project: \url{https://www.worldclim.org/}
#'
#' @examples
#' \donttest{
#' ras <- get_world_pca()
#' # Guarded: the function returns a character string when the data
#' # server is unreachable, so the example also runs offline
#' if (!is.character(ras)) {
#'   print(ras)
#'   terra::plot(ras)
#' }
#' }
get_world_pca <- function() {
  url <- "https://phymaps.nyc3.digitaloceanspaces.com/worldPCA.rds"

  if (!url_exists(url)) {
    return("Can't access data.")
  }

  ras <- readr::read_rds(url)
  terra::unwrap(ras)
}

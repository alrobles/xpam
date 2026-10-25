#' Retrieve Rodentia Range Maps for Mexico
#'
#' Downloads a pre-built \code{\link[terra]{SpatVector}} containing IUCN
#' range polygons for all Rodentia species recorded in Mexico.
#'
#' @return A \code{\link[terra]{SpatVector}} object, or \code{"Can't access data."}
#'   if the remote server is unreachable.
#'
#' @export
#'
#' @source
#' Robles-Fernandez, A. L., Lira-Noriega, A., & Martinez-Meyer, E. (2023).
#' A phylogeny-informed characterisation of global tetrapod traits addresses
#' data gaps and biases.
#' \doi{10.1101/2023.03.04.531098}
#'
#' @examples
#' \donttest{
#' shp <- get_rodentia_mexico()
#' # Guarded: the function returns a character string when the data
#' # server is unreachable, so the example also runs offline
#' if (!is.character(shp)) {
#'   print(shp)
#' }
#' }
get_rodentia_mexico <- function() {
  url <- "https://phymaps.nyc3.digitaloceanspaces.com/rodentia_mexico.rds"

  if (!url_exists(url)) {
    return("Can't access data.")
  }

  shp <- readr::read_rds(url)
  terra::unwrap(shp)
}

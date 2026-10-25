skip_without_gdal <- function() {
  testthat::skip_if_not(has_gdal(), "xpam built without GDAL")
}

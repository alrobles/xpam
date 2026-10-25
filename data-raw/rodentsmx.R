## code to prepare `rodentsmx` dataset goes here
library(mdd)
library(terra)
rodentsmx <- mdd::get_mdd_map(Order = "Rodentia", country = "Mexico", cropCountry = TRUE)
rodentsmx <- terra::wrap(rodentsmx)
readr::write_rds(rodentsmx, "data-raw/rodentsmx.rds", compress = "xz")
usethis::use_data(rodentsmx, overwrite = TRUE, internal = FALSE)

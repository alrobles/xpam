test_that("get_geo_pam rejects non-SpatVector input", {
  expect_error(get_geo_pam(data.frame(x = 1)), "SpatVector")
})

test_that("get_geo_pam rejects invalid resolution", {

  expect_error(get_geo_pam(terra::vect(), res = -1))
  expect_error(get_geo_pam(terra::vect(), res = "a"))
})

test_that("get_geo_pam returns a tibble with expected structure", {
  skip_on_cran()
  skip_if_offline()

  shp <- get_rodentia_mexico()
  skip_if(is.character(shp), "Remote data not available")

  pam <- get_geo_pam(shp, res = 2)

  expect_s3_class(pam, "tbl_df")
  expect_true("x" %in% names(pam))
  expect_true("y" %in% names(pam))
  expect_true(ncol(pam) > 2)

  species_cols <- pam[, -(1:2)]
  expect_true(all(unlist(species_cols) %in% c(0, 1)))
})

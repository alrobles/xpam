test_that("get_world_pca returns SpatRaster or error string", {
  skip_on_cran()
  skip_if_offline()

  result <- get_world_pca()

  if (is.character(result)) {
    expect_equal(result, "Can't access data.")
  } else {
    expect_s4_class(result, "SpatRaster")
    expect_true(terra::nlyr(result) >= 2)
  }
})

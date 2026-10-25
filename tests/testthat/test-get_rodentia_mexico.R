test_that("get_rodentia_mexico returns SpatVector or error string", {
  skip_on_cran()
  skip_if_offline()

  result <- get_rodentia_mexico()

  if (is.character(result)) {
    expect_equal(result, "Can't access data.")
  } else {
    expect_s4_class(result, "SpatVector")
    expect_true("sciname" %in% names(result))
    expect_true(nrow(result) > 0)
  }
})

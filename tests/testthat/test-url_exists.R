test_that("url_exists returns TRUE for reachable URLs", {
  skip_on_cran()
  skip_if_offline()

  expect_true(url_exists("https://cloud.r-project.org"))
})

test_that("url_exists returns FALSE for unreachable URLs", {
  skip_on_cran()

  expect_false(url_exists("https://this-domain-does-not-exist-xyz123.example"))
})

test_that("url_exists prepends https:// when missing", {
  skip_on_cran()
  skip_if_offline()

  expect_true(url_exists("cloud.r-project.org"))
})

test_that("url_exists validates input", {
  expect_error(url_exists(42))
  expect_error(url_exists(NULL))
})

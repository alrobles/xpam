#' Check Whether a URL is Reachable
#'
#' @param url Character string. The URL to check.
#' @param timeout Numeric. Connection timeout in seconds. Default \code{10}.
#'
#' @return \code{TRUE} if the URL responds, \code{FALSE} otherwise.
#' @keywords internal
url_exists <- function(url, timeout = 10) {
  checkmate::assert_string(url)

  if (!grepl("^https?://", url)) {
    url <- paste0("https://", url)
  }

  tryCatch(
    {
      con <- url(url, open = "r")
      on.exit(close(con), add = TRUE)
      TRUE
    },
    error = function(e) FALSE,
    warning = function(w) FALSE
  )
}

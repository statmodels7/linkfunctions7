#' Register the Package's S7 Methods on Load
#'
#' @description
#' Calls `S7::methods_register()`. Without this call, a method that the
#' package registers on a generic of another package does not take effect.
#'
#' @details
#' It matters here for `print` and `plot`: those are S3 generics owned
#' by \pkg{base}, and the `S7::method()` assignments in
#' `methods.R` cannot attach to them until the package is loaded. Without
#' this hook, printing a link object would fall back to the default S7 display.
#'
#' This is the standard R load hook, and it is not called directly.
#'
#' @param ... Ignored; the hook is called by R with the library path and package
#'   name.
#'
#' @return Called for its side effect; the return value is discarded by R.
#'
#' @keywords internal
#' @noRd
.onLoad <- function(...) {
  S7::methods_register()
}

.onLoad <- function(...) {
  # Register S7 methods dynamically when the package is loaded
  S7::methods_register()
}
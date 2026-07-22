#' @title S7 Class for Statistical Link Functions
#'
#' @import S7
#' @description
#' A strictly typed S7 object that encapsulates the metadata of a statistical 
#' link function. The mathematical transformations, including the forward and 
#' inverse functions and their exact analytical derivatives up to the fourth 
#' order, are implemented and registered as S7 generic methods.
#'
#' @details
#' Objects of class \code{link} are instantiated using the S7 object system.
#'
#' The object assumes the following mathematical notation:
#' \itemize{
#'   \item \eqn{\theta}: The response parameter (e.g., probability, mean, dispersion).
#'   \item \eqn{\eta}: The linear predictor (unconstrained scale).
#' }
#' The relationship is defined as \eqn{\eta = g(\theta)} (link function) and
#' \eqn{\theta = g^{-1}(\eta)} (inverse link function).
#'
#' @param link_name A character string identifying the link (e.g., "logit").
#' @param link_bounds A numeric vector of length 2 \code{c(lower, upper)} defining the valid domain for \eqn{\theta}.
#' @param link_params A list or vector of additional parameters required to define the link, or \code{NULL}.
#'
#' @export
link <- S7::new_class(
  name = "link",
  properties = list(
    link_name = S7::class_character,
    link_bounds = S7::class_numeric,
    link_params = S7::class_any
  ),
  
  validator = function(self) {
    # Ensure bounds contain exactly two numeric elements
    if (length(self@link_bounds) != 2) {
      return("Property 'link_bounds' must be a numeric vector of length 2: c(lower, upper).")
    }
    
    # Ensure logical domain definition
    if (self@link_bounds[1] >= self@link_bounds[2]) {
      return("The lower bound must be strictly less than the upper bound.")
    }
  }
)


# A constant of the same length as `v`, missing wherever `v` is missing.
#
# A derivative that reduces to a constant must still report that it does not know
# the answer for an input it was not given. R makes this easy to get wrong:
# `NA^0` is 1, so `theta^(lambda - 2)` silently turns a missing parameter into a
# number as soon as `lambda` is 2.
const_like <- function(v, value) {
  out <- rep(value, length(v))
  out[is.na(v)] <- NA_real_
  out
}

# Carry missingness from the input over to a computed result, for the same reason
# as const_like(): an expression whose exponent happens to vanish stops depending
# on its argument, and loses the missingness along with it.
na_from <- function(r, v) {
  r[is.na(v)] <- NA_real_
  r
}

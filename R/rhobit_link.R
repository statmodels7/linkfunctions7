RhobitLink <- S7::new_class(
  name = "RhobitLink",
  parent = link
)

# --- Methods for RhobitLink ---

# Forward and inverse link functions
S7::method(linkfun, RhobitLink) <- function(x, theta) atanh(theta)
S7::method(linkinv, RhobitLink) <- function(x, eta) tanh(eta)

# Exact analytical derivatives of the link function (wrt theta)
S7::method(dlinkfun, RhobitLink) <- function(x, theta) {
  1 / (1 - theta^2)
}
S7::method(d2linkfun, RhobitLink) <- function(x, theta) {
  (2 * theta) / ((1 - theta^2)^2)
}
S7::method(d3linkfun, RhobitLink) <- function(x, theta) {
  (2 + 6 * (theta^2)) / ((1 - theta^2)^3)
}
S7::method(d4linkfun, RhobitLink) <- function(x, theta) {
  (24 * theta * (1 + theta^2)) / ((1 - theta^2)^4)
}

# Exact analytical derivatives of the inverse link function (wrt eta)
# We elegantly evaluate all these as polynomials of t = tanh(eta)
S7::method(dlinkinv, RhobitLink) <- function(x, eta) {
  t <- tanh(eta); 1 - t^2
}
S7::method(d2linkinv, RhobitLink) <- function(x, eta) {
  t <- tanh(eta); -2 * t * (1 - t^2)
}
S7::method(d3linkinv, RhobitLink) <- function(x, eta) {
  t <- tanh(eta); -2 + 8 * (t^2) - 6 * (t^4)
}
S7::method(d4linkinv, RhobitLink) <- function(x, eta) {
  t <- tanh(eta); 16 * t - 40 * (t^3) + 24 * (t^5)
}

#' @title The Rhobit (Fisher's z) Link Function
#'
#' @include generics.R
#' @include link_class.R
#' @description
#' Creates an S7 object of class \code{link} implementing the Rhobit transformation, 
#' also known as Fisher's z-transformation. This link function rigorously maps the 
#' open interval \code{c(-1, 1)} to the real line \code{c(-Inf, Inf)}. It is primarily 
#' used for modeling correlation coefficients or other bounded parameters that are 
#' symmetrically constrained.
#'
#' @details
#' The Rhobit link is defined mathematically using the inverse hyperbolic tangent function:
#' \eqn{\eta = \text{arctanh}(\theta) = \frac{1}{2} \log\left(\frac{1 + \theta}{1 - \theta}\right)}.
#'
#' The inverse link is the hyperbolic tangent function:
#' \eqn{\theta = \tanh(\eta) = \frac{\exp(2\eta) - 1}{\exp(2\eta) + 1}}.
#'
#' The valid mathematical domain of \eqn{\theta} is exactly \code{c(-1, 1)}.
#'
#' @return An S7 object of class \code{RhobitLink} (inheriting from \code{link}) containing the transformation functions
#' and their exact analytical derivatives up to the fourth order.
#'
#' @seealso \code{\link{link}}, \code{\link{logit_link}}
#' @export
rhobit_link <- function() {
  RhobitLink(
    link_name = "rhobit",
    link_bounds = c(-1, 1),
    link_params = NULL
  )
}
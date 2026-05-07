SqrtLink <- S7::new_class(
  name = "SqrtLink",
  parent = link
)

# --- Methods for SqrtLink ---

# Forward and inverse link functions
S7::method(linkfun, SqrtLink) <- function(x, theta) sqrt(theta)
S7::method(linkinv, SqrtLink) <- function(x, eta) eta^2

# Exact analytical derivatives of the link function (wrt theta)
S7::method(dlinkfun, SqrtLink) <- function(x, theta) {
  1 / (2 * sqrt(theta))
}
S7::method(d2linkfun, SqrtLink) <- function(x, theta) {
  -1 / (4 * (theta^1.5))
}
S7::method(d3linkfun, SqrtLink) <- function(x, theta) {
  3 / (8 * (theta^2.5))
}
S7::method(d4linkfun, SqrtLink) <- function(x, theta) {
  -15 / (16 * (theta^3.5))
}

# Exact analytical derivatives of the inverse link function (wrt eta)
# 3rd and 4th derivatives uniquely vanish to exactly 0 for this quadratic form.
S7::method(dlinkinv, SqrtLink) <- function(x, eta) 2 * eta
S7::method(d2linkinv, SqrtLink) <- function(x, eta) rep(2, length(eta))
S7::method(d3linkinv, SqrtLink) <- function(x, eta) rep(0, length(eta))
S7::method(d4linkinv, SqrtLink) <- function(x, eta) rep(0, length(eta))

#' @title The Square Root Link Function
#'
#' @include generics.R
#' @include link_class.R
#' @description
#' Creates an S7 object of class \code{link} implementing the Square Root transformation.
#' This is a specific case of the Power link family (with \eqn{\lambda = 0.5}) and is 
#' prominently used for modeling count data (e.g., Poisson regression) as a 
#' variance-stabilizing transformation.
#'
#' @details
#' The Square Root link is mathematically defined as \eqn{\eta = \sqrt{\theta}}.
#' Consequently, the inverse link is derived as \eqn{\theta = \eta^2}.
#'
#' Unlike the Log link, this transformation allows \eqn{\theta} to reach 0 exactly.
#' While the inverse function (\eqn{\eta^2}) is mathematically valid for negative 
#' values of \eqn{\eta}, in the specific context of this link function, the linear 
#' predictor \eqn{\eta} is typically constrained to be non-negative. This restriction 
#' preserves a strictly one-to-one mapping with \eqn{\theta}.
#'
#' The strict mathematical domain for \eqn{\theta} is \code{c(0, Inf)}.
#'
#' @return An S7 object of class \code{SqrtLink} (inheriting from \code{link}) containing the transformation functions
#' and their exact analytical derivatives up to the fourth order.
#'
#' @seealso \code{\link{link}}, \code{\link{power_link}}, \code{\link{log_link}}
#' @export
sqrt_link <- function() {
  SqrtLink(
    link_name = "sqrt",
    link_bounds = c(0, Inf),
    
    # The standard square root link requires no additional parameters
    link_params = NULL
  )
}
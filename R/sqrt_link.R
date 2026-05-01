#' @title The Square Root Link Function
#'
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
#' @return An S7 object of class \code{link} containing the transformation functions
#' and their exact analytical derivatives up to the fourth order.
#'
#' @seealso \code{\link{link}}, \code{\link{power_link}}, \code{\link{log_link}}
#'
#' @export
sqrt_link <- function() {
  link(
    link_name = "sqrt",
    link_bounds = c(0, Inf),
    
    # The standard square root link requires no additional parameters
    link_params = NULL,
    
    # Forward and inverse link functions
    linkfun = function(theta) sqrt(theta),
    linkinv = function(eta) eta^2,
    
    # Exact analytical derivatives of the link function (wrt theta)
    dlinkfun  = function(theta) {
      1 / (2 * sqrt(theta))
    },
    d2linkfun = function(theta) {
      -1 / (4 * (theta^1.5))
    },
    d3linkfun = function(theta) {
      3 / (8 * (theta^2.5))
    },
    d4linkfun = function(theta) {
      -15 / (16 * (theta^3.5))
    },
    
    # Exact analytical derivatives of the inverse link function (wrt eta)
    # 3rd and 4th derivatives uniquely vanish to exactly 0 for this quadratic form.
    dlinkinv  = function(eta) 2 * eta,
    d2linkinv = function(eta) rep(2, length(eta)),
    d3linkinv = function(eta) rep(0, length(eta)),
    d4linkinv = function(eta) rep(0, length(eta))
  )
}
#' @title The Inverse Square Link Function
#'
#' @description
#' Creates an S7 object of class \code{link} implementing the Inverse Square transformation.
#' This is the canonical link function for the Inverse Gaussian distribution.
#'
#' @details
#' The Inverse Square link is defined mathematically as \eqn{\eta = 1 / \theta^2}.
#' Consequently, the inverse link function is derived as \eqn{\theta = 1 / \sqrt{\eta}}.
#'
#' This specific link function is predominantly utilized in Generalized Linear Models (GLMs)
#' assuming an Inverse Gaussian response distribution. In such frameworks, the variance
#' is proportional to the cube of the mean (\eqn{\text{Var}(Y) \propto \mu^3}).
#'
#' \strong{Domain and Optimization Constraints:}
#' Both the parameter \eqn{\theta} and the linear predictor \eqn{\eta} must be strictly
#' positive. The valid mathematical domain for \eqn{\theta} is \code{c(0, Inf)}. During
#' optimization routines (e.g., Fisher Scoring or Newton-Raphson), extreme care must be
#' taken to ensure the linear predictor \eqn{\eta > 0}. Evaluating the inverse link or
#' its derivatives at non-positive values of \eqn{\eta} will inevitably result in \code{NaN}s
#' due to fractional powers and square root operations.
#'
#' @return An S7 object of class \code{link} containing the transformation functions
#' and their exact analytical derivatives up to the fourth order.
#'
#' @seealso \code{\link{link}}, \code{\link{inverse_link}}
#'
#' @export
inverse_sq_link <- function() {
  link(
    link_name = "inverse_sq",
    link_bounds = c(0, Inf),
    
    # The inverse square link requires no additional mathematical parameters
    link_params = NULL,
    
    # Forward and inverse link functions
    linkfun = function(theta) 1 / (theta^2),
    linkinv = function(eta) 1 / sqrt(eta),
    
    # Exact analytical derivatives of the link function (wrt theta)
    dlinkfun  = function(theta) -2 / (theta^3),
    d2linkfun = function(theta)  6 / (theta^4),
    d3linkfun = function(theta) -24 / (theta^5),
    d4linkfun = function(theta) 120 / (theta^6),
    
    # Exact analytical derivatives of the inverse link function (wrt eta)
    # Utilizing explicit fractions to optimize numeric evaluations
    dlinkinv  = function(eta) -1 / (2 * eta^1.5),
    d2linkinv = function(eta)  3 / (4 * eta^2.5),
    d3linkinv = function(eta) -15 / (8 * eta^3.5),
    d4linkinv = function(eta) 105 / (16 * eta^4.5)
  )
}
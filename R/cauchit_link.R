#' @title The Cauchit Link Function
#'
#' @description
#' Creates an S7 object of class \code{link} implementing the Cauchit transformation.
#' This link function rigorously maps the open interval \code{c(0, 1)} to the real line 
#' by utilizing the quantile function of the standard Cauchy distribution.
#'
#' @details
#' The Cauchit link is defined mathematically as \eqn{\eta = \tan(\pi(\theta - 0.5))}, 
#' which corresponds perfectly to \code{qcauchy(theta)}.
#' The inverse link is the standard Cauchy CDF \eqn{\theta = \frac{1}{\pi} \arctan(\eta) + 0.5},
#' computed via \code{pcauchy(eta)}.
#'
#' \strong{Heavy Tails:} Unlike the Logit or Probit links, the Cauchit link has 
#' exceedingly heavier tails. This makes it particularly robust and useful for modeling 
#' binary data where the probability approaches 0 or 1 very slowly, or when the dataset 
#' contains severe outliers that might disproportionately influence the fit of 
#' light-tailed link functions.
#'
#' The strictly valid mathematical domain for \eqn{\theta} is \code{c(0, 1)}.
#'
#' @return An S7 object of class \code{link} containing the transformation functions
#' and their exact analytical derivatives up to the fourth order.
#'
#' @seealso \code{\link{link}}, \code{\link{logit_link}}, \code{\link{probit_link}}
#' @importFrom stats qcauchy pcauchy
#' @export
cauchit_link <- function() {
  link(
    link_name = "cauchit",
    link_bounds = c(0, 1),
    
    link_params = NULL,
    
    # Forward and inverse link functions utilizing native C-level implementation
    linkfun = function(theta) stats::qcauchy(theta),
    linkinv = function(eta) stats::pcauchy(eta),
    
    # Exact analytical derivatives of the link function (wrt theta)
    # We compute eta locally to maximize computational efficiency 
    # and express higher-order derivatives as elegant polynomials of eta.
    dlinkfun  = function(theta) {
      eta <- stats::qcauchy(theta)
      pi * (1 + eta^2)
    },
    d2linkfun = function(theta) {
      eta <- stats::qcauchy(theta)
      2 * (pi^2) * eta * (1 + eta^2)
    },
    d3linkfun = function(theta) {
      eta <- stats::qcauchy(theta)
      2 * (pi^3) * (1 + eta^2) * (1 + 3 * (eta^2))
    },
    d4linkfun = function(theta) {
      eta <- stats::qcauchy(theta)
      8 * (pi^4) * eta * (1 + eta^2) * (2 + 3 * (eta^2))
    },
    
    # Exact analytical derivatives of the inverse link function (wrt eta)
    # Derived purely from the Cauchy probability density function
    dlinkinv  = function(eta) {
      1 / (pi * (1 + eta^2))
    },
    d2linkinv = function(eta) {
      -2 * eta / (pi * ((1 + eta^2)^2))
    },
    d3linkinv = function(eta) {
      2 * (3 * (eta^2) - 1) / (pi * ((1 + eta^2)^3))
    },
    d4linkinv = function(eta) {
      24 * eta * (1 - eta^2) / (pi * ((1 + eta^2)^4))
    }
  )
}
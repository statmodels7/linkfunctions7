#' @title The Probit Link Function
#'
#' @description
#' Creates an S7 object of class \code{link} implementing the Probit transformation.
#' This link function rigorously relies on the cumulative distribution function (CDF) 
#' of the standard Normal distribution. It is widely used in Generalized Linear Models 
#' (GLMs) for binary data, often justified by a latent normal variable interpretation.
#'
#' @details
#' The Probit link is mathematically defined as \eqn{\eta = \Phi^{-1}(\theta)}, where 
#' \eqn{\Phi^{-1}} is the quantile function of the standard normal distribution (\code{qnorm}).
#' The inverse link is \eqn{\theta = \Phi(\eta)}, the standard normal CDF (\code{pnorm}).
#'
#' Similarly to the \code{logit} link, the Probit is symmetric around \eqn{\theta = 0.5} 
#' (where \eqn{\eta = 0}). However, the tails of the Normal distribution approach 0 
#' and 1 faster than the Logistic distribution.
#'
#' The strictly mathematical domain of \eqn{\theta} is \code{c(0, 1)}.
#'
#' @return An S7 object of class \code{link} containing the transformation functions
#' and their exact analytical derivatives up to the fourth order.
#'
#' @seealso \code{\link{link}}, \code{\link{logit_link}}, \code{\link{cauchit_link}}
#' @importFrom stats qnorm pnorm dnorm
#' @export
probit_link <- function() {
  link(
    link_name = "probit",
    link_bounds = c(0, 1),
    
    link_params = NULL,
    
    linkfun = function(theta) stats::qnorm(theta),
    linkinv = function(eta) stats::pnorm(eta),
    
    # Exact analytical derivatives of the link function (wrt theta)
    # We calculate eta = qnorm(theta) and phi = dnorm(eta) locally 
    # to significantly optimize computational operations.
    dlinkfun  = function(theta) {
      eta <- stats::qnorm(theta)
      1 / stats::dnorm(eta)
    },
    d2linkfun = function(theta) {
      eta <- stats::qnorm(theta)
      phi <- stats::dnorm(eta)
      eta / (phi^2)
    },
    d3linkfun = function(theta) {
      eta <- stats::qnorm(theta)
      phi <- stats::dnorm(eta)
      (1 + 2 * (eta^2)) / (phi^3)
    },
    d4linkfun = function(theta) {
      eta <- stats::qnorm(theta)
      phi <- stats::dnorm(eta)
      (7 * eta + 6 * (eta^3)) / (phi^4)
    },
    
    # Exact analytical derivatives of the inverse link function (wrt eta)
    # They rely exclusively on standard normal density properties.
    dlinkinv  = function(eta) {
      stats::dnorm(eta)
    },
    d2linkinv = function(eta) {
      phi <- stats::dnorm(eta)
      -eta * phi
    },
    d3linkinv = function(eta) {
      phi <- stats::dnorm(eta)
      (eta^2 - 1) * phi
    },
    d4linkinv = function(eta) {
      phi <- stats::dnorm(eta)
      (3 * eta - eta^3) * phi
    }
  )
}
#' @title The Logit Link Function
#'
#' @description
#' Creates an S7 object of class \code{link} implementing the Logit (log-odds) transformation.
#' This is the canonical link function for the success probability parameter of the Bernoulli 
#' and Binomial distributions and serves as the foundation of Logistic Regression.
#'
#' @details
#' The Logit link is defined mathematically as \eqn{\eta = \log(\frac{\theta}{1 - \theta})}.
#' The inverse link is the standard logistic function (sigmoid): 
#' \eqn{\theta = \frac{1}{1 + \exp(-\eta)}}.
#'
#' The link is perfectly symmetric around \eqn{\theta = 0.5} (which corresponds to \eqn{\eta = 0}). 
#' It elegantly interprets the linear predictor as the log-odds of the event probability.
#'
#' The mathematical domain of \eqn{\theta} is \code{c(0, 1)}.
#'
#' \strong{Implementation Details:} 
#' This function internally delegates to R's highly optimized native functions 
#' \code{stats::qlogis} and \code{stats::plogis}. This ensures maximum numerical stability 
#' and precision, especially when evaluating probabilities exceedingly close to the 
#' boundaries of 0 and 1.
#'
#' @return An S7 object of class \code{link} containing the transformation functions
#' and their exact analytical derivatives up to the fourth order.
#'
#' @seealso \code{\link{link}}, \code{\link{probit_link}}, \code{\link{cloglog_link}}
#' @importFrom stats qlogis plogis
#' @export
logit_link <- function() {
  link(
    link_name = "logit",
    link_bounds = c(0, 1),
    
    # The logit link requires no additional mathematical parameters
    link_params = NULL,
    
    # Forward and inverse link functions using native R C-level functions
    linkfun = function(theta) stats::qlogis(theta),
    linkinv = function(eta) stats::plogis(eta),
    
    # Exact analytical derivatives of the link function (wrt theta)
    dlinkfun  = function(theta) 1 / (theta * (1 - theta)),
    d2linkfun = function(theta) (2 * theta - 1) / ((theta * (1 - theta))^2),
    d3linkfun = function(theta) 2 / (theta^3) + 2 / ((1 - theta)^3),
    d4linkfun = function(theta) -6 / (theta^4) + 6 / ((1 - theta)^4),
    
    # Exact analytical derivatives of the inverse link function (wrt eta)
    # Expressed elegantly as polynomials of the probability p to maximize performance
    dlinkinv  = function(eta) {
      p <- stats::plogis(eta)
      p * (1 - p)
    },
    d2linkinv = function(eta) {
      p <- stats::plogis(eta)
      p * (1 - p) * (1 - 2 * p)
    },
    d3linkinv = function(eta) {
      p <- stats::plogis(eta)
      p * (1 - p) * (1 - 6 * p + 6 * p^2)
    },
    d4linkinv = function(eta) {
      p <- stats::plogis(eta)
      p * (1 - p) * (1 - 14 * p + 36 * p^2 - 24 * p^3)
    }
  )
}
#' @title The Log-Log Link Function
#'
#' @description
#' Creates an S7 object of class \code{link} implementing the Log-Log transformation.
#' This link is the asymmetric opposite of the \code{cloglog} link and is often used 
#' in survival analysis (e.g., Gompertz models) or for modeling binary data with 
#' specific asymmetry requirements.
#'
#' @details
#' The Log-Log link is mathematically defined as \eqn{\eta = -\log(-\log(\theta))}.
#' Consequently, the inverse link is derived as \eqn{\theta = \exp(-\exp(-\eta))}.
#'
#' \strong{Asymmetry:} Unlike the symmetric Logit or Probit links, the Log-Log link 
#' is highly asymmetric. It is particularly suitable for modeling events where the 
#' probability approaches 0 very slowly but approaches 1 sharply. This represents 
#' the mathematical reverse of the \code{cloglog} link.
#'
#' The strict mathematical domain for \eqn{\theta} is \code{c(0, 1)}.
#'
#' @return An S7 object of class \code{link} containing the transformation functions
#' and their exact analytical derivatives up to the fourth order.
#'
#' @seealso \code{\link{link}}, \code{\link{cloglog_link}}, \code{\link{logit_link}}
#'
#' @export
loglog_link <- function() {
  link(
    link_name = "loglog",
    link_bounds = c(0, 1),
    
    link_params = NULL,
    
    linkfun = function(theta) -log(-log(theta)),
    linkinv = function(eta) exp(-exp(-eta)),
    
    # Exact analytical derivatives of the link function (wrt theta)
    dlinkfun  = function(theta) {
      -1 / (theta * log(theta))
    },
    d2linkfun = function(theta) {
      l <- log(theta)
      (1 + l) / (theta^2 * (l^2))
    },
    d3linkfun = function(theta) {
      l <- log(theta)
      -(2 + 3 * l + 2 * (l^2)) / (theta^3 * (l^3))
    },
    d4linkfun = function(theta) {
      l <- log(theta)
      (6 + 12 * l + 11 * (l^2) + 6 * (l^3)) / (theta^4 * (l^4))
    },
    
    # Exact analytical derivatives of the inverse link function (wrt eta)
    # Utilizing the term z = exp(-eta) to evaluate derivatives as polynomials,
    # thus maximizing computational performance.
    dlinkinv  = function(eta) { z <- exp(-eta); exp(-z) * z },
    d2linkinv = function(eta) { z <- exp(-eta); exp(-z) * (z^2 - z) },
    d3linkinv = function(eta) { z <- exp(-eta); exp(-z) * (z^3 - 3 * z^2 + z) },
    d4linkinv = function(eta) { z <- exp(-eta); exp(-z) * (z^4 - 6 * z^3 + 7 * z^2 - z) }
  )
}
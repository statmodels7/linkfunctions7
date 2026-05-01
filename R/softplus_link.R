#' @title The Softplus Link Function
#'
#' @description
#' Creates an S7 object of class \code{link} implementing the Softplus transformation.
#' The Softplus function is a smooth approximation of the rectifier (ReLU) and ensures 
#' the parameter \eqn{\theta} remains positive. Unlike the Log link, which implies an 
#' exponential relationship throughout, the Softplus link asymptotically approaches 
#' linearity for large positive values of the linear predictor.
#'
#' @param a A numeric value specifying the scaling parameter (smoothness/steepness). 
#' Must be strictly positive. Defaults to 1.
#'
#' @details
#' The Softplus link describes the relationship where the response parameter \eqn{\theta} 
#' is the Softplus of the linear predictor \eqn{\eta}.
#'
#' Mathematically:
#' \itemize{
#'   \item Inverse Link (Softplus): \eqn{\theta = \frac{1}{a} \log(1 + \exp(a \eta))}
#'   \item Link Function: \eqn{\eta = \frac{1}{a} \log(\exp(a \theta) - 1)}
#' }
#'
#' \strong{Behavior:}
#' For large negative \eqn{\eta}, \eqn{\theta \approx 0}.
#' For large positive \eqn{\eta}, \eqn{\theta \approx \eta} (linear behavior), whereas 
#' a Log link would imply \eqn{\theta = \exp(\eta)} (exponential behavior).
#'
#' \strong{Numerical Stability:}
#' The inverse link implementation intelligently utilizes conditional algebraic logic 
#' to ensure robust numerical stability for large positive values of \eqn{\eta}, 
#' entirely avoiding precision overflow.
#'
#' The mathematical domain of \eqn{\theta} is \code{c(0, Inf)}.
#'
#' @return An S7 object of class \code{link} containing the transformation functions,
#' their exact analytical derivatives up to the fourth order, and the parameter \code{a}.
#'
#' @seealso \code{\link{link}}, \code{\link{log_link}}, \code{\link{identity_link}}
#' @importFrom stats plogis
#' @export
softplus_link <- function(a = 1) {
  if (a <= 0) {
    stop("Scale parameter 'a' must be strictly greater than 0.")
  }

  link(
    link_name = paste0("softplus(a=", round(a, 5), ")"),
    link_bounds = c(0, Inf),
    link_params = list(a = a),
    
    # Forward and inverse link functions with numerical stability guards
    linkfun = function(theta) {
      log(expm1(a * theta)) / a
    },
    linkinv = function(eta) {
      # "Log-Sum-Exp" algebraic trick to avoid computing exp() of large numbers
      pmax(0, eta) + log(1 + exp(-abs(a * eta))) / a
    },
    
    # Exact analytical derivatives of the link function (wrt theta)
    # We use x = expm1(a * theta) to express derivatives safely and efficiently
    dlinkfun  = function(theta) {
      x <- expm1(a * theta); (x + 1) / x
    },
    d2linkfun = function(theta) {
      x <- expm1(a * theta); -a * (x + 1) / (x^2)
    },
    d3linkfun = function(theta) {
      x <- expm1(a * theta); (a^2) * (x + 1) * (x + 2) / (x^3)
    },
    d4linkfun = function(theta) {
      x <- expm1(a * theta); -(a^3) * (x + 1) * (x^2 + 6 * x + 6) / (x^4)
    },
    
    # Exact analytical derivatives of the inverse link function (wrt eta)
    # They cleanly resolve into polynomials of the logistic function p
    dlinkinv  = function(eta) {
      stats::plogis(a * eta)
    },
    d2linkinv = function(eta) {
      p <- stats::plogis(a * eta); a * p * (1 - p)
    },
    d3linkinv = function(eta) {
      p <- stats::plogis(a * eta); (a^2) * p * (1 - p) * (1 - 2 * p)
    },
    d4linkinv = function(eta) {
      p <- stats::plogis(a * eta); (a^3) * p * (1 - p) * (1 - 6 * p + 6 * (p^2))
    }
  )
}
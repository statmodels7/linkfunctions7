#' @title The General Bounded Link Function
#'
#' @description
#' Creates an S7 object of class \code{link} that maps a constrained interval to the 
#' whole real line. By specifying \code{lwr} and \code{upr}, this function dynamically 
#' constructs a doubly bounded (scaled logit), lower bounded (shifted log), upper bounded, 
#' or unbounded (identity) link function.
#'
#' @param lwr Numeric or \code{NULL}. The lower bound of the interval.
#' @param upr Numeric or \code{NULL}. The upper bound of the interval.
#'
#' @details
#' \strong{Doubly Bounded (\code{lwr} and \code{upr} provided):}
#' Transforms \eqn{\theta} by normalizing it to \code{c(0, 1)} via 
#' \eqn{p = \frac{\theta - \text{lwr}}{\text{upr} - \text{lwr}}}, and then applying the logit function.
#' 
#' \strong{Lower Bounded (\code{lwr} provided, \code{upr = NULL}):}
#' Defined as \eqn{\eta = \log(\theta - \text{lwr})}, with inverse \eqn{\theta = \exp(\eta) + \text{lwr}}.
#' 
#' \strong{Upper Bounded (\code{lwr = NULL}, \code{upr} provided):}
#' Defined as \eqn{\eta = \log(\text{upr} - \theta)}, with inverse \eqn{\theta = \text{upr} - \exp(\eta)}.
#' 
#' \strong{Unbounded (\code{lwr = NULL}, \code{upr = NULL}):}
#' Returns the standard \code{\link{identity_link}}.
#'
#' @return An S7 object of class \code{link} containing the transformation functions
#' and their exact analytical derivatives up to the fourth order.
#'
#' @seealso \code{\link{link}}, \code{\link{logit_link}}, \code{\link{log_link}}
#' @importFrom stats qlogis plogis
#' @export
bounded_link <- function(lwr = NULL, upr = NULL) {
  
  # Case 0: Unbounded link effectively reduces to Identity
  if (is.null(lwr) && is.null(upr)) {
    return(identity_link())
  }
  
  safe_exp <- function(eta) {
    pmax(exp(eta), .Machine$double.eps)
  }
  
  # Case 1: Doubly Bounded
  if (!is.null(lwr) && !is.null(upr)) {
    if (lwr >= upr) stop("Lower bound 'lwr' must be strictly less than upper bound 'upr'.")
    
    W <- upr - lwr
    
    return(link(
      link_name = paste0("bounded(lwr=", lwr, ", upr=", upr, ")"),
      link_bounds = c(lwr, upr),
      link_params = list(lwr = lwr, upr = upr),
      
      linkfun = function(theta) {
        p <- (theta - lwr) / W
        stats::qlogis(p)
      },
      linkinv = function(eta) {
        lwr + W * stats::plogis(eta)
      },
      
      dlinkfun  = function(theta) {
        p <- (theta - lwr) / W
        1 / (W * p * (1 - p))
      },
      d2linkfun = function(theta) {
        p <- (theta - lwr) / W
        (2 * p - 1) / ((W^2) * (p^2) * ((1 - p)^2))
      },
      d3linkfun = function(theta) {
        p <- (theta - lwr) / W
        (2 / (p^3) + 2 / ((1 - p)^3)) / (W^3)
      },
      d4linkfun = function(theta) {
        p <- (theta - lwr) / W
        (-6 / (p^4) + 6 / ((1 - p)^4)) / (W^4)
      },
      
      dlinkinv  = function(eta) { p <- stats::plogis(eta); W * p * (1 - p) },
      d2linkinv = function(eta) { p <- stats::plogis(eta); W * p * (1 - p) * (1 - 2 * p) },
      d3linkinv = function(eta) { p <- stats::plogis(eta); W * p * (1 - p) * (1 - 6 * p + 6 * (p^2)) },
      d4linkinv = function(eta) { p <- stats::plogis(eta); W * p * (1 - p) * (1 - 14 * p + 36 * (p^2) - 24 * (p^3)) }
    ))
  }
  
  # Case 2: Lower Bounded
  if (!is.null(lwr) && is.null(upr)) {
    return(link(
      link_name = paste0("lower_bounded(lwr=", lwr, ")"),
      link_bounds = c(lwr, Inf),
      link_params = list(lwr = lwr),
      
      linkfun = function(theta) log(theta - lwr),
      linkinv = function(eta) lwr + safe_exp(eta),
      
      dlinkfun  = function(theta)  1 / (theta - lwr),
      d2linkfun = function(theta) -1 / ((theta - lwr)^2),
      d3linkfun = function(theta)  2 / ((theta - lwr)^3),
      d4linkfun = function(theta) -6 / ((theta - lwr)^4),
      
      dlinkinv  = safe_exp,
      d2linkinv = safe_exp,
      d3linkinv = safe_exp,
      d4linkinv = safe_exp
    ))
  }
  
  # Case 3: Upper Bounded
  if (is.null(lwr) && !is.null(upr)) {
    return(link(
      link_name = paste0("upper_bounded(upr=", upr, ")"),
      link_bounds = c(-Inf, upr),
      link_params = list(upr = upr),
      
      linkfun = function(theta) log(upr - theta),
      linkinv = function(eta) upr - safe_exp(eta),
      
      dlinkfun  = function(theta) {
        -1 / (upr - theta)
      },
      d2linkfun = function(theta) {
        -1 / ((upr - theta)^2)
      },
      d3linkfun = function(theta) {
        -2 / ((upr - theta)^3)
      },
      d4linkfun = function(theta) {
        -6 / ((upr - theta)^4)
      },
      
      dlinkinv  = function(eta) -safe_exp(eta),
      d2linkinv = function(eta) -safe_exp(eta),
      d3linkinv = function(eta) -safe_exp(eta),
      d4linkinv = function(eta) -safe_exp(eta)
    ))
  }
}
#' @title The Logarithmic Link Function
#'
#' @description
#' Creates an S7 object of class \code{link} implementing the natural logarithm transformation.
#' This is the canonical link function for the mean of the Poisson distribution and is 
#' widely used for modeling count data or non-negative continuous data with multiplicative effects.
#'
#' @details
#' The Log link is defined mathematically as \eqn{\eta = \log(\theta)}.
#' The inverse link is the exponential function \eqn{\theta = \exp(\eta)}.
#'
#' A remarkable mathematical property of this link is that the inverse function is its 
#' own derivative. Therefore, the parameter \eqn{\theta} and all its derivatives with 
#' respect to \eqn{\eta} are equal to \eqn{\exp(\eta)}.
#'
#' The valid mathematical domain of \eqn{\theta} is \code{c(0, Inf)}. 
#'
#' \strong{Numerical Stability:}
#' During the evaluation of the inverse link and its derivatives, the result is bounded 
#' from below by \code{.Machine$double.eps}. This prevents numerical underflow to 
#' exactly zero when \eqn{\eta} is a large negative number, which would otherwise 
#' produce \code{Inf} when subsequently calculating \eqn{1/\theta}.
#'
#' @return An S7 object of class \code{link} containing the transformation functions
#' and their exact analytical derivatives up to the fourth order.
#'
#' @seealso \code{\link{link}}, \code{\link{inverse_link}}
#'
#' @export
log_link <- function() {
  
  # Auxiliary function to ensure numerical stability against underflow to zero
  safe_exp <- function(eta) {
    pmax(exp(eta), .Machine$double.eps)
  }
  
  link(
    link_name = "log",
    link_bounds = c(0, Inf),
    
    # The log link requires no additional mathematical parameters
    link_params = NULL,
    
    # Forward and inverse link functions
    linkfun = function(theta) log(theta),
    linkinv = safe_exp,
    
    # Exact analytical derivatives of the link function (wrt theta)
    dlinkfun  = function(theta)  1 / theta,
    d2linkfun = function(theta) -1 / (theta^2),
    d3linkfun = function(theta)  2 / (theta^3),
    d4linkfun = function(theta) -6 / (theta^4),
    
    # Exact analytical derivatives of the inverse link function (wrt eta)
    # d^k/deta^k exp(eta) = exp(eta) for all k > 0
    dlinkinv  = safe_exp,
    d2linkinv = safe_exp,
    d3linkinv = safe_exp,
    d4linkinv = safe_exp
  )
}
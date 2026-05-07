LogLink <- S7::new_class(
  name = "LogLink",
  parent = link
)

# --- Methods for LogLink ---

# Forward and inverse link functions
S7::method(linkfun, LogLink) <- function(x, theta) log(theta)
S7::method(linkinv, LogLink) <- function(x, eta) pmax(exp(eta), .Machine$double.eps)

# Exact analytical derivatives of the link function (wrt theta)
S7::method(dlinkfun, LogLink) <- function(x, theta)  1 / theta
S7::method(d2linkfun, LogLink) <- function(x, theta) -1 / (theta^2)
S7::method(d3linkfun, LogLink) <- function(x, theta)  2 / (theta^3)
S7::method(d4linkfun, LogLink) <- function(x, theta) -6 / (theta^4)

# Exact analytical derivatives of the inverse link function (wrt eta)
# d^k/deta^k exp(eta) = exp(eta) for all k > 0, bounded for numerical stability
S7::method(dlinkinv, LogLink) <- function(x, eta) pmax(exp(eta), .Machine$double.eps)
S7::method(d2linkinv, LogLink) <- function(x, eta) pmax(exp(eta), .Machine$double.eps)
S7::method(d3linkinv, LogLink) <- function(x, eta) pmax(exp(eta), .Machine$double.eps)
S7::method(d4linkinv, LogLink) <- function(x, eta) pmax(exp(eta), .Machine$double.eps)

#' @title The Logarithmic Link Function
#'
#' @include generics.R
#' @include link_class.R
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
#' @return An S7 object of class \code{LogLink} (inheriting from \code{link}) containing the transformation functions
#' and their exact analytical derivatives up to the fourth order.
#'
#' @seealso \code{\link{link}}, \code{\link{inverse_link}}
#' @export
log_link <- function() {
  LogLink(
    link_name = "log",
    link_bounds = c(0, Inf),
    
    # The log link requires no additional mathematical parameters
    link_params = NULL
  )
}
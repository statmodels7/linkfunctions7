#' @title S7 Class for the Logarithmic Link
#'
#' @description
#' The class \code{\link{log_link}} instantiates.
#'
#' @return An S7 object of class \code{LogLink}, inheriting from
#'   \code{\link{link}}.
#'
#' @seealso \code{\link{log_link}}, the constructor users call.
#' @keywords internal
LogLink <- S7::new_class(
  name = "LogLink",
  parent = link
)

# --- Methods for LogLink ---

# Forward and inverse link functions
S7::method(linkfun, LogLink) <- function(x, theta) log(theta)
S7::method(linkinv, LogLink) <- function(x, eta) exp_floored(eta)

# Exact analytical derivatives of the link function (wrt theta)
S7::method(dlinkfun, LogLink) <- function(x, theta)  1 / theta
S7::method(d2linkfun, LogLink) <- function(x, theta) -1 / (theta^2)
S7::method(d3linkfun, LogLink) <- function(x, theta)  2 / (theta^3)
S7::method(d4linkfun, LogLink) <- function(x, theta) -6 / (theta^4)

# Exact analytical derivatives of the inverse link function (wrt eta)
# d^k/deta^k exp(eta) = exp(eta) for all k > 0, floored for numerical stability
S7::method(dlinkinv, LogLink) <- function(x, eta) exp_floored(eta)
S7::method(d2linkinv, LogLink) <- function(x, eta) exp_floored(eta)
S7::method(d3linkinv, LogLink) <- function(x, eta) exp_floored(eta)
S7::method(d4linkinv, LogLink) <- function(x, eta) exp_floored(eta)

#' @title The Logarithmic Link Function
#'
#' @include generics.R
#' @include link_class.R
#' @description
#' The log link \eqn{\eta = \log\theta} on \eqn{(0, \infty)}, with
#' inverse \eqn{\theta = e^\eta}; the canonical link for a positive
#' parameter.
#' @details
#' The Log link is defined mathematically as \eqn{\eta = \log(\theta)}.
#' The inverse link is the exponential function \eqn{\theta = \exp(\eta)}.
#'
#' The inverse function of this link is its 
#' own derivative. Therefore, the parameter \eqn{\theta} and all its derivatives with 
#' respect to \eqn{\eta} are equal to \eqn{\exp(\eta)}.
#'
#' The valid mathematical domain of \eqn{\theta} is \code{c(0, Inf)}. 
#'
#' \strong{Numerical Stability:}
#' The inverse link and its derivatives are bounded below by \code{exp_floor},
#' which is \code{.Machine$double.xmin^0.25}, about \code{1.2e-77}. This prevents
#' underflow to exactly zero for large negative \eqn{\eta}, which would produce
#' \code{Inf} when the forward derivatives divide by \eqn{\theta}; the fourth of
#' them divides by \eqn{\theta^4}, and that is what sets the value. The floor is
#' low enough that \eqn{\theta} is exact down to \eqn{\eta \approx -177}.
#'
#' @return An S7 object of class \code{LogLink} (inheriting from \code{link}) containing the transformation functions
#' and their exact analytical derivatives up to the fourth order.
#'
#' @examples
#' lk <- log_link()
#' lk
#'
#' theta <- c(0.5, 1, 10)
#' eta <- linkfun(lk, theta)
#' eta
#' linkinv(lk, eta)
#'
#' # the exponential is its own derivative, so every inverse derivative agrees
#' dlinkinv(lk, eta)
#' d4linkinv(lk, eta)
#'
#' # forward derivatives to fourth order
#' linkderiv(lk, theta, order = 4)
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
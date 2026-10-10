#' @title S7 Class for the Logarithmic Link
#'
#' @description
#' Carries the log transformation \eqn{\eta = \log\theta} on \eqn{(0, \infty)},
#' with inverse \eqn{\theta = e^{\eta}}. It is the standard link for a
#' positive parameter such as a scale or a rate.
#'
#' The inverse is floored at [exp_floor()], so a parameter reported by this link
#' is never exactly zero and can be divided into.
#'
#' @param link_name A character string naming the link, set by the
#'   constructor and shown by `print()`.
#' @param link_bounds A length-two numeric vector, the open interval in
#'   which the parameter lies. Set by the constructor; see Value for this link's.
#' @param link_params A list of the link's own parameters, empty where it has
#'   none. Set by the constructor.
#'
#' @return An S7 object of class `LogLink`, inheriting from [link()] and
#'   carrying its three properties `link_name`, `link_bounds` and
#'   `link_params`. Its `link_bounds` are `c(0, Inf)` and it carries no link parameters.
#'
#' @section Methods:
#' Twelve methods are registered on this class: [linkfun()] and [linkinv()],
#' and the five derivative orders in each direction, [dlinkfun()] through
#' [d5linkfun()] going out and [dlinkinv()] through [d5linkinv()] coming
#' back. The forward derivatives are \eqn{(-1)^{k-1}(k-1)!\,\theta^{-k}},
#' written out. Every inverse derivative is [exp_floored()] of \eqn{\eta},
#' the exponential being its own derivative to every order, and the floor is
#' what keeps \eqn{\theta} strictly inside \eqn{(0, \infty)}.
#'
#' @aliases linkfun.LogLink
#' @aliases linkinv.LogLink
#' @aliases dlinkfun.LogLink
#' @aliases d2linkfun.LogLink
#' @aliases d3linkfun.LogLink
#' @aliases d4linkfun.LogLink
#' @aliases d5linkfun.LogLink
#' @aliases dlinkinv.LogLink
#' @aliases d2linkinv.LogLink
#' @aliases d3linkinv.LogLink
#' @aliases d4linkinv.LogLink
#' @aliases d5linkinv.LogLink
#'
#' @seealso [log_link()], the constructor users call.
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
S7::method(d5linkfun, LogLink) <- function(x, theta) 24 / (theta^5)

# Exact analytical derivatives of the inverse link function (wrt eta)
# d^k/deta^k exp(eta) = exp(eta) for all k > 0, floored for numerical stability
S7::method(dlinkinv, LogLink) <- function(x, eta) exp_floored(eta)
S7::method(d2linkinv, LogLink) <- function(x, eta) exp_floored(eta)
S7::method(d3linkinv, LogLink) <- function(x, eta) exp_floored(eta)
S7::method(d4linkinv, LogLink) <- function(x, eta) exp_floored(eta)
S7::method(d5linkinv, LogLink) <- function(x, eta) exp_floored(eta)

#' @title The Logarithmic Link Function
#'
#' @include generics.R
#' @include link_class.R
#' @description
#' The log link \eqn{\eta = \log\theta} on \eqn{(0, \infty)}, with
#' inverse \eqn{\theta = e^\eta}; the canonical link for a positive
#' parameter.
#' @details
#' The log link is \eqn{\eta = \log(\theta)}, and its inverse is the
#' exponential function \eqn{\theta = \exp(\eta)}.
#'
#' The exponential is its own derivative, so every derivative of the inverse
#' link with respect to \eqn{\eta} equals \eqn{\exp(\eta)}.
#'
#' The domain of \eqn{\theta} is \eqn{(0, \infty)}.
#'
#' The inverse link and its derivatives are bounded below by [exp_floor()],
#' \eqn{(24/x_{\max})^{1/4} \approx 1.9 \times 10^{-77}}, where \eqn{x_{\max}}
#' is the largest finite double. The floor prevents an underflow to exactly
#' zero for a large negative \eqn{\eta}, which would make the forward
#' derivatives infinite; its value is set by the fourth of them, which divides
#' by \eqn{\theta^4}. The floor is low enough that \eqn{\theta} is exact down
#' to \eqn{\eta \approx -176.7}.
#'
#' @return An S7 object of class `LogLink`, inheriting from [link()], whose
#'   methods compute the link, its inverse and their derivatives to the fifth
#'   order.
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
#' @seealso [link()], [inverse_link()]
#' @export
log_link <- function() {
  LogLink(
    link_name = "log",
    link_bounds = c(0, Inf),
    
    # The log link requires no additional mathematical parameters
    link_params = NULL
  )
}

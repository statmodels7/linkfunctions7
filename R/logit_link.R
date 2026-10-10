#' @title S7 Class for the Logit Link
#'
#' @description
#' Carries the logit transformation \eqn{\eta = \log(\theta/(1-\theta))} on
#' \eqn{(0, 1)}, with inverse \eqn{\theta = 1/(1+e^{-\eta})}. It is the
#' canonical link for a probability, and its linear predictor is the log-odds.
#'
#' Every inverse derivative is a polynomial in \eqn{\theta} itself, the first
#' being \eqn{\theta(1-\theta)}, the variance of a Bernoulli.
#'
#' @param link_name A character string naming the link, set by the
#'   constructor and shown by `print()`.
#' @param link_bounds A length-two numeric vector, the open interval in
#'   which the parameter lies. Set by the constructor; see Value for this link's.
#' @param link_params A list of the link's own parameters, empty where it has
#'   none. Set by the constructor.
#'
#' @return An S7 object of class `LogitLink`, inheriting from [link()] and
#'   carrying its three properties `link_name`, `link_bounds` and
#'   `link_params`. Its `link_bounds` are `c(0, 1)` and it carries no link parameters.
#'
#' @section Methods:
#' Twelve methods are registered on this class: [linkfun()] and [linkinv()],
#' and the five derivative orders in each direction, [dlinkfun()] through
#' [d5linkfun()] going out and [dlinkinv()] through [d5linkinv()] coming
#' back. `linkfun()` and `linkinv()` delegate to `stats::qlogis()` and
#' `stats::plogis()`, which remain accurate near both boundaries. The ten
#' derivatives come from a compiled kernel, one call per order and
#' direction. The inverse derivatives are evaluated in \eqn{\theta},
#' \eqn{1 - \theta} and their product, each computed from \eqn{\eta}, so
#' that they keep their relative accuracy in both tails.
#'
#' @aliases linkfun.LogitLink
#' @aliases linkinv.LogitLink
#' @aliases dlinkfun.LogitLink
#' @aliases d2linkfun.LogitLink
#' @aliases d3linkfun.LogitLink
#' @aliases d4linkfun.LogitLink
#' @aliases d5linkfun.LogitLink
#' @aliases dlinkinv.LogitLink
#' @aliases d2linkinv.LogitLink
#' @aliases d3linkinv.LogitLink
#' @aliases d4linkinv.LogitLink
#' @aliases d5linkinv.LogitLink
#'
#' @seealso [logit_link()], the constructor users call.
#' @keywords internal
LogitLink <- S7::new_class(
  name = "LogitLink",
  parent = link
)

# --- Methods for LogitLink ---

# Forward and inverse link functions using native R C-level functions
S7::method(linkfun, LogitLink) <- function(x, theta) stats::qlogis(theta)
S7::method(linkinv, LogitLink) <- function(x, eta) stats::plogis(eta)

# Exact analytical derivatives of the link function (wrt theta)
S7::method(dlinkfun, LogitLink) <- function(x, theta) lk_logit_fwd_cpp(theta, 1L)
S7::method(d2linkfun, LogitLink) <- function(x, theta) lk_logit_fwd_cpp(theta, 2L)
S7::method(d3linkfun, LogitLink) <- function(x, theta) lk_logit_fwd_cpp(theta, 3L)
S7::method(d4linkfun, LogitLink) <- function(x, theta) lk_logit_fwd_cpp(theta, 4L)
S7::method(d5linkfun, LogitLink) <- function(x, theta) lk_logit_fwd_cpp(theta, 5L)

# Exact analytical derivatives of the inverse link function (wrt eta).
#
# Polynomials in the probability p itself, evaluated by logistic_poly() inside
# the kernel. They are shared with the doubly bounded link, which scales them by
# the interval width, and with the softplus, which uses them one order down.
S7::method(dlinkinv, LogitLink) <- function(x, eta) lk_logit_inv_cpp(eta, 1L)
S7::method(d2linkinv, LogitLink) <- function(x, eta) lk_logit_inv_cpp(eta, 2L)
S7::method(d3linkinv, LogitLink) <- function(x, eta) lk_logit_inv_cpp(eta, 3L)
S7::method(d4linkinv, LogitLink) <- function(x, eta) lk_logit_inv_cpp(eta, 4L)
S7::method(d5linkinv, LogitLink) <- function(x, eta) lk_logit_inv_cpp(eta, 5L)

#' @title The Logit Link Function
#'
#' @include generics.R
#' @include link_class.R
#' @description
#' The logit link \eqn{\eta = \log(\theta/(1-\theta))} on \eqn{(0, 1)},
#' with inverse \eqn{\theta = 1/(1+e^{-\eta})}; the canonical link for a
#' probability, whose linear predictor is the log-odds.
#' @details
#' The logit link is \eqn{\eta = \log(\theta/(1 - \theta))}, and its inverse
#' is the logistic function \eqn{\theta = 1/(1 + \exp(-\eta))}.
#'
#' The link is symmetric about \eqn{\theta = 1/2}, where \eqn{\eta = 0}, and
#' the linear predictor is the log-odds of the event probability. The domain
#' of \eqn{\theta} is \eqn{(0, 1)}.
#'
#' The implementation calls `stats::qlogis()` and `stats::plogis()`, which
#' remain accurate near both boundaries.
#'
#' @return An S7 object of class `LogitLink`, inheriting from [link()], whose
#'   methods compute the link, its inverse and their derivatives to the fifth
#'   order.
#'
#' @examples
#' lk <- logit_link()
#' lk
#'
#' p <- c(0.1, 0.5, 0.9)
#' eta <- linkfun(lk, p)      # log-odds
#' eta
#' linkinv(lk, eta)
#'
#' # the inverse derivatives are polynomials in p: the first is the variance
#' # of a Bernoulli, p(1 - p)
#' dlinkinv(lk, 0)
#'
#' # all four orders at once
#' vapply(1:4, function(k) linkinvderiv(lk, 0, order = k), numeric(1))
#'
#' # every mathematical property is checkable
#' check_link(lk)
#'
#' @seealso [link()], [probit_link()], [cloglog_link()]
#' @importFrom stats qlogis plogis
#' @export
logit_link <- function() {
  LogitLink(
    link_name = "logit",
    link_bounds = c(0, 1),

    # The logit link requires no additional mathematical parameters
    link_params = NULL
  )
}

#' @title S7 Class for the LogLog Link
#'
#' @description
#' Carries the log-log transformation \eqn{\eta = -\log(-\log\theta)} on
#' \eqn{(0, 1)}, with inverse \eqn{\theta = \exp(-e^{-\eta})}.
#'
#' It is the mirror image of [cloglog_link()] about \eqn{\theta = 1/2}:
#' reflecting one link's parameter gives the other's, so it approaches zero
#' faster than one.
#'
#' @param link_name A character string naming the link, set by the
#'   constructor and shown by `print()`.
#' @param link_bounds A length-two numeric vector, the open interval in
#'   which the parameter lies. Set by the constructor; see Value for this link's.
#' @param link_params A list of the link's own parameters, empty where it has
#'   none. Set by the constructor.
#'
#' @return An S7 object of class `LogLogLink`, inheriting from [link()] and
#'   carrying its three properties `link_name`, `link_bounds` and
#'   `link_params`. Its `link_bounds` are `c(0, 1)` and it carries no link parameters.
#'
#' @section Methods:
#' Twelve methods are registered on this class: [linkfun()] and [linkinv()],
#' and the five derivative orders in each direction, [dlinkfun()] through
#' [d5linkfun()] going out and [dlinkinv()] through [d5linkinv()] coming
#' back. `linkfun()` and `linkinv()` are the elementary
#' \eqn{-\log(-\log\theta)} and \eqn{\exp(-\exp(-\eta))}; the ten
#' derivatives come from a compiled kernel.
#'
#' @aliases linkfun.LogLogLink
#' @aliases linkinv.LogLogLink
#' @aliases dlinkfun.LogLogLink
#' @aliases d2linkfun.LogLogLink
#' @aliases d3linkfun.LogLogLink
#' @aliases d4linkfun.LogLogLink
#' @aliases d5linkfun.LogLogLink
#' @aliases dlinkinv.LogLogLink
#' @aliases d2linkinv.LogLogLink
#' @aliases d3linkinv.LogLogLink
#' @aliases d4linkinv.LogLogLink
#' @aliases d5linkinv.LogLogLink
#'
#' @seealso [loglog_link()], the constructor users call.
#' @keywords internal
LogLogLink <- S7::new_class(
  name = "LogLogLink",
  parent = link
)

# --- Methods for LogLogLink ---

S7::method(linkfun, LogLogLink) <- function(x, theta) -log(-log(theta))
S7::method(linkinv, LogLogLink) <- function(x, eta) exp(-exp(-eta))

# Exact analytical derivatives of the link function (wrt theta)
S7::method(dlinkfun, LogLogLink) <- function(x, theta) lk_loglog_fwd_cpp(theta, 1L)
S7::method(d2linkfun, LogLogLink) <- function(x, theta) lk_loglog_fwd_cpp(theta, 2L)

S7::method(d3linkfun, LogLogLink) <- function(x, theta) lk_loglog_fwd_cpp(theta, 3L)
S7::method(d4linkfun, LogLogLink) <- function(x, theta) lk_loglog_fwd_cpp(theta, 4L)
S7::method(d5linkfun, LogLogLink) <- function(x, theta) lk_loglog_fwd_cpp(theta, 5L)

# Exact analytical derivatives of the inverse link function (wrt eta)
# Utilizing the term z = exp(-eta) to evaluate derivatives as polynomials,
# thus maximizing computational performance.
S7::method(dlinkinv, LogLogLink) <- function(x, eta) lk_loglog_inv_cpp(eta, 1L)
S7::method(d2linkinv, LogLogLink) <- function(x, eta) lk_loglog_inv_cpp(eta, 2L)
S7::method(d3linkinv, LogLogLink) <- function(x, eta) lk_loglog_inv_cpp(eta, 3L)
S7::method(d4linkinv, LogLogLink) <- function(x, eta) lk_loglog_inv_cpp(eta, 4L)
S7::method(d5linkinv, LogLogLink) <- function(x, eta) lk_loglog_inv_cpp(eta, 5L)

#' @title The Log-Log Link Function
#'
#' @include generics.R
#' @include link_class.R
#' @description
#' The log-log link \eqn{\eta = -\log(-\log\theta)} on \eqn{(0, 1)},
#' with inverse \eqn{\theta = \exp(-e^{-\eta})}; the mirror image of
#' [cloglog_link()].
#' @details
#' The log-log link is \eqn{\eta = -\log(-\log\theta)}, with inverse
#' \eqn{\theta = \exp(-\exp(-\eta))}.
#'
#' Unlike the logit and the probit, the link is asymmetric: \eqn{\theta}
#' approaches 0 faster than 1, which is the mirror image of [cloglog_link()].
#' The domain of \eqn{\theta} is \eqn{(0, 1)}.
#'
#' @return An S7 object of class `LogLogLink`, inheriting from [link()], whose
#'   methods compute the link, its inverse and their derivatives to the fifth
#'   order.
#'
#' @examples
#' lk <- loglog_link()
#' lk
#'
#' p <- c(0.1, 0.5, 0.9)
#' eta <- linkfun(lk, p)
#' eta
#' linkinv(lk, eta)
#'
#' # the mirror image of the cloglog link
#' linkinv(loglog_link(), 1)
#' 1 - linkinv(cloglog_link(), -1)
#'
#' linkderiv(lk, 0.5, order = 3)
#'
#' @seealso [link()], [cloglog_link()], [logit_link()]
#' @export
loglog_link <- function() {
  LogLogLink(
    link_name = "loglog",
    link_bounds = c(0, 1),
    link_params = NULL
  )
}

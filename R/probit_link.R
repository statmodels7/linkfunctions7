#' @title S7 Class for the Probit Link
#'
#' @description
#' Carries the probit transformation \eqn{\eta = \Phi^{-1}(\theta)} on
#' \eqn{(0, 1)}, with \eqn{\Phi} the standard normal distribution function and
#' the inverse \eqn{\theta = \Phi(\eta)}.
#'
#' It is symmetric about \eqn{\theta = 1/2}, like the logit, and it reaches
#' its bounds faster, because the tails of the normal distribution are lighter
#' than those of the logistic.
#'
#' @param link_name A character string naming the link, set by the
#'   constructor and shown by `print()`.
#' @param link_bounds A length-two numeric vector, the open interval in
#'   which the parameter lies. Set by the constructor; see Value for this link's.
#' @param link_params A list of the link's own parameters, empty where it has
#'   none. Set by the constructor.
#'
#' @return An S7 object of class `ProbitLink`, inheriting from [link()] and
#'   carrying its three properties `link_name`, `link_bounds` and
#'   `link_params`. Its `link_bounds` are `c(0, 1)` and it carries no link parameters.
#'
#' @section Methods:
#' Twelve methods are registered on this class: [linkfun()] and [linkinv()],
#' and the five derivative orders in each direction, [dlinkfun()] through
#' [d5linkfun()] going out and [dlinkinv()] through [d5linkinv()] coming
#' back. `linkfun()` and `linkinv()` delegate to `stats::qnorm()` and
#' `stats::pnorm()`. The ten derivatives come from a compiled kernel, one
#' call per order and direction.
#'
#' @aliases linkfun.ProbitLink
#' @aliases linkinv.ProbitLink
#' @aliases dlinkfun.ProbitLink
#' @aliases d2linkfun.ProbitLink
#' @aliases d3linkfun.ProbitLink
#' @aliases d4linkfun.ProbitLink
#' @aliases d5linkfun.ProbitLink
#' @aliases dlinkinv.ProbitLink
#' @aliases d2linkinv.ProbitLink
#' @aliases d3linkinv.ProbitLink
#' @aliases d4linkinv.ProbitLink
#' @aliases d5linkinv.ProbitLink
#'
#' @seealso [probit_link()], the constructor users call.
#' @keywords internal
ProbitLink <- S7::new_class(
  name = "ProbitLink",
  parent = link
)

# --- Methods for ProbitLink ---

S7::method(linkfun, ProbitLink) <- function(x, theta) stats::qnorm(theta)
S7::method(linkinv, ProbitLink) <- function(x, eta) stats::pnorm(eta)

# Exact analytical derivatives of the link function (wrt theta)
# We calculate eta = qnorm(theta) and phi = dnorm(eta) locally 
# to significantly optimize computational operations.
S7::method(dlinkfun, ProbitLink) <- function(x, theta) lk_probit_fwd_cpp(theta, 1L)
S7::method(d2linkfun, ProbitLink) <- function(x, theta) lk_probit_fwd_cpp(theta, 2L)
S7::method(d3linkfun, ProbitLink) <- function(x, theta) lk_probit_fwd_cpp(theta, 3L)
S7::method(d4linkfun, ProbitLink) <- function(x, theta) lk_probit_fwd_cpp(theta, 4L)
S7::method(d5linkfun, ProbitLink) <- function(x, theta) lk_probit_fwd_cpp(theta, 5L)

# Exact analytical derivatives of the inverse link function (wrt eta)
# They rely exclusively on standard normal density properties.
S7::method(dlinkinv, ProbitLink) <- function(x, eta) lk_probit_inv_cpp(eta, 1L)
S7::method(d2linkinv, ProbitLink) <- function(x, eta) lk_probit_inv_cpp(eta, 2L)
S7::method(d3linkinv, ProbitLink) <- function(x, eta) lk_probit_inv_cpp(eta, 3L)
S7::method(d4linkinv, ProbitLink) <- function(x, eta) lk_probit_inv_cpp(eta, 4L)
S7::method(d5linkinv, ProbitLink) <- function(x, eta) lk_probit_inv_cpp(eta, 5L)

#' @title The Probit Link Function
#'
#' @include generics.R
#' @include link_class.R
#' @description
#' The probit link \eqn{\eta = \Phi^{-1}(\theta)} on \eqn{(0, 1)}, with
#' \eqn{\Phi} the standard normal distribution function.
#' @details
#' The probit link is \eqn{\eta = \Phi^{-1}(\theta)}, where \eqn{\Phi^{-1}}
#' is the quantile function of the standard normal distribution (`qnorm()`).
#' The inverse link is the standard normal distribution function,
#' \eqn{\theta = \Phi(\eta)} (`pnorm()`).
#'
#' Like the logit, the probit link is symmetric about \eqn{\theta = 1/2},
#' where \eqn{\eta = 0}. The tails of the normal distribution are lighter than
#' those of the logistic, so \eqn{\theta} approaches 0 and 1 faster than under
#' the logit link.
#'
#' The domain of \eqn{\theta} is \eqn{(0, 1)}.
#'
#' @return An S7 object of class `ProbitLink`, inheriting from [link()], whose
#'   methods compute the link, its inverse and their derivatives to the fifth
#'   order.
#'
#' @examples
#' lk <- probit_link()
#' lk
#'
#' p <- c(0.1, 0.5, 0.9)
#' eta <- linkfun(lk, p)      # standard normal quantiles
#' eta
#' linkinv(lk, eta)
#'
#' # the first inverse derivative is the standard normal density
#' dlinkinv(lk, 0)
#' dnorm(0)
#'
#' # probit tails approach 0 and 1 faster than logit ones
#' linkinv(probit_link(), 3)
#' linkinv(logit_link(), 3)
#'
#' @seealso [link()], [logit_link()], [cauchit_link()]
#' @importFrom stats qnorm pnorm dnorm
#' @export
probit_link <- function() {
  ProbitLink(
    link_name = "probit",
    link_bounds = c(0, 1),
    link_params = NULL
  )
}

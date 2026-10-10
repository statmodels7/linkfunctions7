#' @title S7 Class for the Cauchit Link
#'
#' @description
#' Carries the cauchit transformation \eqn{\eta = \tan(\pi(\theta - 1/2))} on
#' \eqn{(0, 1)}, the Cauchy quantile function, with inverse
#' \eqn{\theta = 1/2 + \arctan(\eta)/\pi}.
#'
#' Its tails are heavier than those of the logit and the probit, so an extreme
#' linear predictor moves the probability less, and a few observations are
#' less able to drive the fit to a boundary.
#'
#' @param link_name A character string naming the link, set by the
#'   constructor and shown by `print()`.
#' @param link_bounds A length-two numeric vector, the open interval in
#'   which the parameter lies. Set by the constructor; see Value for this link's.
#' @param link_params A list of the link's own parameters, empty where it has
#'   none. Set by the constructor.
#'
#' @return An S7 object of class `CauchitLink`, inheriting from [link()] and
#'   carrying its three properties `link_name`, `link_bounds` and
#'   `link_params`. Its `link_bounds` are `c(0, 1)` and it carries no link parameters.
#'
#' @section Methods:
#' Twelve methods are registered on this class: [linkfun()] and [linkinv()],
#' and the five derivative orders in each direction, [dlinkfun()] through
#' [d5linkfun()] going out and [dlinkinv()] through [d5linkinv()] coming
#' back. `linkfun()` and `linkinv()` delegate to `stats::qcauchy()` and
#' `stats::pcauchy()`, which stay accurate in both tails. The ten
#' derivatives come from a compiled kernel, one call per order and
#' direction.
#'
#' @aliases linkfun.CauchitLink
#' @aliases linkinv.CauchitLink
#' @aliases dlinkfun.CauchitLink
#' @aliases d2linkfun.CauchitLink
#' @aliases d3linkfun.CauchitLink
#' @aliases d4linkfun.CauchitLink
#' @aliases d5linkfun.CauchitLink
#' @aliases dlinkinv.CauchitLink
#' @aliases d2linkinv.CauchitLink
#' @aliases d3linkinv.CauchitLink
#' @aliases d4linkinv.CauchitLink
#' @aliases d5linkinv.CauchitLink
#'
#' @seealso [cauchit_link()], the constructor users call.
#' @keywords internal
CauchitLink <- S7::new_class(
  name = "CauchitLink",
  parent = link
)

# --- Methods for CauchitLink ---

# Forward and inverse link functions utilizing native C-level implementation
S7::method(linkfun, CauchitLink) <- function(x, theta) stats::qcauchy(theta)
S7::method(linkinv, CauchitLink) <- function(x, eta) stats::pcauchy(eta)

# Exact analytical derivatives of the link function (wrt theta)
# We compute eta locally to maximize computational efficiency 
# and express higher-order derivatives as elegant polynomials of eta.
S7::method(dlinkfun, CauchitLink) <- function(x, theta) lk_cauchit_fwd_cpp(theta, 1L)
S7::method(d2linkfun, CauchitLink) <- function(x, theta) lk_cauchit_fwd_cpp(theta, 2L)
S7::method(d3linkfun, CauchitLink) <- function(x, theta) lk_cauchit_fwd_cpp(theta, 3L)
S7::method(d4linkfun, CauchitLink) <- function(x, theta) lk_cauchit_fwd_cpp(theta, 4L)
S7::method(d5linkfun, CauchitLink) <- function(x, theta) lk_cauchit_fwd_cpp(theta, 5L)

# Exact analytical derivatives of the inverse link function (wrt eta)
# Derived purely from the Cauchy probability density function
S7::method(dlinkinv, CauchitLink) <- function(x, eta) lk_cauchit_inv_cpp(eta, 1L)
S7::method(d2linkinv, CauchitLink) <- function(x, eta) lk_cauchit_inv_cpp(eta, 2L)
S7::method(d3linkinv, CauchitLink) <- function(x, eta) lk_cauchit_inv_cpp(eta, 3L)
S7::method(d4linkinv, CauchitLink) <- function(x, eta) lk_cauchit_inv_cpp(eta, 4L)
S7::method(d5linkinv, CauchitLink) <- function(x, eta) lk_cauchit_inv_cpp(eta, 5L)

#' @title The Cauchit Link Function
#'
#' @include generics.R
#' @include link_class.R
#' @description
#' The cauchit link \eqn{\eta = \tan(\pi(\theta - 1/2))} on
#' \eqn{(0, 1)}, the Cauchy quantile function; heavier-tailed than the
#' logit or the probit.
#' @details
#' The cauchit link is \eqn{\eta = \tan(\pi(\theta - 1/2))}, computed as
#' `qcauchy(theta)`. The inverse link is the standard Cauchy distribution
#' function \eqn{\theta = \arctan(\eta)/\pi + 1/2}, computed as
#' `pcauchy(eta)`.
#'
#' The Cauchy distribution has heavier tails than the logistic and the normal,
#' so the probability approaches 0 and 1 more slowly than under the logit or
#' the probit link. The link is therefore less sensitive to observations whose
#' linear predictor is extreme.
#'
#' The domain of \eqn{\theta} is \eqn{(0, 1)}.
#'
#' @return An S7 object of class `CauchitLink`, inheriting from [link()], whose
#'   methods compute the link, its inverse and their derivatives to the fifth
#'   order.
#'
#' @examples
#' lk <- cauchit_link()
#' lk
#'
#' p <- c(0.1, 0.5, 0.9)
#' eta <- linkfun(lk, p)
#' eta
#' linkinv(lk, eta)
#'
#' # heavy tails: the same eta is far less extreme than under a logit
#' linkinv(cauchit_link(), 5)
#' linkinv(logit_link(), 5)
#'
#' dlinkinv(lk, 0)            # 1 / pi
#'
#' @seealso [link()], [logit_link()], [probit_link()]
#' @importFrom stats qcauchy pcauchy
#' @export
cauchit_link <- function() {
  CauchitLink(
    link_name = "cauchit",
    link_bounds = c(0, 1),
    link_params = NULL
  )
}

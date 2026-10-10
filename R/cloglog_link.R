#' @title S7 Class for the Complementary Log-Log Link
#'
#' @description
#' Carries the complementary log-log transformation
#' \eqn{\eta = \log(-\log(1-\theta))} on \eqn{(0, 1)}, with inverse
#' \eqn{\theta = 1 - \exp(-e^{\eta})}.
#'
#' Unlike the logit and the probit it is asymmetric about \eqn{\theta = 1/2},
#' approaching one faster than zero, so it is the link of a
#' proportional-hazards model for a binary outcome.
#'
#' @param link_name A character string naming the link, set by the
#'   constructor and shown by `print()`.
#' @param link_bounds A length-two numeric vector, the open interval in
#'   which the parameter lies. Set by the constructor; see Value for this link's.
#' @param link_params A list of the link's own parameters, empty where it has
#'   none. Set by the constructor.
#'
#' @return An S7 object of class `ClogLogLink`, inheriting from [link()] and
#'   carrying its three properties `link_name`, `link_bounds` and
#'   `link_params`. Its `link_bounds` are `c(0, 1)` and it carries no link parameters.
#'
#' @section Methods:
#' Twelve methods are registered on this class: [linkfun()] and [linkinv()],
#' and the five derivative orders in each direction, [dlinkfun()] through
#' [d5linkfun()] going out and [dlinkinv()] through [d5linkinv()] coming
#' back. `linkfun()` is written through `log1p()` rather than as
#' \eqn{\log(-\log(1 - \theta))}, which rounds to \eqn{-\infty} for a small
#' \eqn{\theta} where the true value is finite and representable. The ten
#' derivatives come from a compiled kernel.
#'
#' @aliases linkfun.ClogLogLink
#' @aliases linkinv.ClogLogLink
#' @aliases dlinkfun.ClogLogLink
#' @aliases d2linkfun.ClogLogLink
#' @aliases d3linkfun.ClogLogLink
#' @aliases d4linkfun.ClogLogLink
#' @aliases d5linkfun.ClogLogLink
#' @aliases dlinkinv.ClogLogLink
#' @aliases d2linkinv.ClogLogLink
#' @aliases d3linkinv.ClogLogLink
#' @aliases d4linkinv.ClogLogLink
#' @aliases d5linkinv.ClogLogLink
#'
#' @seealso [cloglog_link()], the constructor users call.
#' @keywords internal
ClogLogLink <- S7::new_class(
  name = "ClogLogLink",
  parent = link
)

# --- Methods for ClogLogLink ---

S7::method(linkfun, ClogLogLink) <- function(x, theta) {
  # log1p(-theta), not log(1 - theta). For small theta the subtraction rounds to
  # exactly 1, its logarithm to exactly 0, and the answer to -Inf -- while the
  # true value, log(-log(1-theta)) ~ log(theta), is perfectly representable:
  # at theta = 1.9e-77 it is -176.66. linkinv() reaches that far down on purpose
  # (see its floor below), so linkfun has to come back from there.
  log(-log1p(-theta))
}

S7::method(linkinv, ClogLogLink) <- function(x, eta) {
  # Kept strictly inside (0, 1). The lower floor is exp_floor rather than
  # .Machine$double.eps because for very negative eta the true value is
  # -expm1(-exp(eta)) ~ exp(eta), which is representable far below eps; flooring
  # at eps threw away every value below eta = -36. -expm1() also keeps the small
  # end accurate, where 1 - exp(-z) cancels.
  pmax(pmin(-expm1(-exp(eta)), 1 - .Machine$double.eps), exp_floor)
}

# Exact analytical derivatives of the link function (wrt theta)
S7::method(dlinkfun, ClogLogLink) <- function(x, theta) lk_cloglog_fwd_cpp(theta, 1L)
S7::method(d2linkfun, ClogLogLink) <- function(x, theta) lk_cloglog_fwd_cpp(theta, 2L)
S7::method(d3linkfun, ClogLogLink) <- function(x, theta) lk_cloglog_fwd_cpp(theta, 3L)
S7::method(d4linkfun, ClogLogLink) <- function(x, theta) lk_cloglog_fwd_cpp(theta, 4L)
S7::method(d5linkfun, ClogLogLink) <- function(x, theta) lk_cloglog_fwd_cpp(theta, 5L)

# Exact analytical derivatives of the inverse link function (wrt eta).
#
# Written as T_k(eta) = exp(k*eta - exp(eta)), which sidesteps the "Inf * 0" NaN
# that comes of forming exp(k*eta) and exp(-exp(eta)) separately. The floor that
# used to be applied to z here was a no-op -- subtracting 2.2e-16 from eta cannot
# change exp(eta) -- so it is gone.
S7::method(dlinkinv, ClogLogLink) <- function(x, eta) lk_cloglog_inv_cpp(eta, 1L)
S7::method(d2linkinv, ClogLogLink) <- function(x, eta) lk_cloglog_inv_cpp(eta, 2L)
S7::method(d3linkinv, ClogLogLink) <- function(x, eta) lk_cloglog_inv_cpp(eta, 3L)
S7::method(d4linkinv, ClogLogLink) <- function(x, eta) lk_cloglog_inv_cpp(eta, 4L)
S7::method(d5linkinv, ClogLogLink) <- function(x, eta) lk_cloglog_inv_cpp(eta, 5L)

#' @title The Complementary Log-Log (ClogLog) Link Function
#'
#' @include generics.R
#' @include link_class.R
#' @description
#' The complementary log-log link \eqn{\eta = \log(-\log(1-\theta))} on
#' \eqn{(0, 1)}, with inverse \eqn{\theta = 1 - \exp(-e^\eta)};
#' asymmetric about \eqn{\theta = 1/2}.
#' @details
#' The complementary log-log link is \eqn{\eta = \log(-\log(1 - \theta))},
#' with inverse \eqn{\theta = 1 - \exp(-\exp(\eta))}.
#'
#' Unlike the logit and the probit, the link is asymmetric: \eqn{\theta}
#' approaches 1 faster than 0. The inverse link is the distribution function
#' of the Gumbel distribution for minima. The link is used in discrete-time
#' survival analysis, where it gives a proportional-hazards model, and for
#' rare events.
#'
#' The domain of \eqn{\theta} is \eqn{(0, 1)}.
#'
#' @return An S7 object of class `ClogLogLink`, inheriting from [link()], whose
#'   methods compute the link, its inverse and their derivatives to the fifth
#'   order.
#'
#' @examples
#' lk <- cloglog_link()
#' lk
#'
#' p <- c(0.1, 0.5, 0.9)
#' eta <- linkfun(lk, p)
#' eta
#' linkinv(lk, eta)
#'
#' # asymmetric: it reaches 1 sharply and 0 slowly, the mirror of loglog
#' linkinv(cloglog_link(), c(-2, 2))
#' linkinv(loglog_link(),  c(-2, 2))
#'
#' d2linkinv(lk, 0)
#'
#' @seealso [link()], [logit_link()], [loglog_link()]
#' @export
cloglog_link <- function() {
  ClogLogLink(
    link_name = "cloglog",
    link_bounds = c(0, 1),
    link_params = NULL
  )
}

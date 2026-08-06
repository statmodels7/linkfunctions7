#' @title S7 Class for the Complementary Log-Log Link
#'
#' @description
#' The class \code{\link{cloglog_link}} instantiates.
#'
#' @return An S7 object of class \code{ClogLogLink}, inheriting from
#'   \code{\link{link}}.
#'
#' @seealso \code{\link{cloglog_link}}, the constructor users call.
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

#' @title The Complementary Log-Log (ClogLog) Link Function
#'
#' @include generics.R
#' @include link_class.R
#' @description
#' The complementary log-log link \eqn{\eta = \log(-\log(1-\theta))} on
#' \eqn{(0, 1)}, with inverse \eqn{\theta = 1 - \exp(-e^\eta)};
#' asymmetric about \eqn{\theta = 1/2}.
#' @details
#' The ClogLog link is defined mathematically as \eqn{\eta = \log(-\log(1 - \theta))}.
#' Consequently, the inverse link is derived as \eqn{\theta = 1 - \exp(-\exp(\eta))}.
#'
#' Unlike the symmetric Logit and Probit links, the ClogLog link lacks symmetry. 
#' It is fundamentally related to the Extreme Value (Gumbel) distribution and is 
#' frequently utilized in discrete-time survival analysis (proportional hazards models) 
#' as well as for modeling rare events.
#'
#' The strictly valid mathematical domain for \eqn{\theta} is \code{c(0, 1)}.
#'
#' @return An S7 object of class \code{ClogLogLink} (inheriting from \code{link})
#' containing the transformation functions and their exact analytical derivatives
#' up to the fourth order.
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
#' # asymmetric: it reaches 1 slowly and 0 sharply, the mirror of loglog
#' linkinv(cloglog_link(), c(-2, 2))
#' linkinv(loglog_link(),  c(-2, 2))
#'
#' d2linkinv(lk, 0)
#'
#' @seealso \code{\link{link}}, \code{\link{logit_link}}, \code{\link{loglog_link}}
#' @export
cloglog_link <- function() {
  ClogLogLink(
    link_name = "cloglog",
    link_bounds = c(0, 1),
    link_params = NULL
  )
}

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
  log(-log(1 - theta))
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
S7::method(dlinkfun, ClogLogLink) <- function(x, theta) {
  val <- 1 - theta
  L <- log(val)
  -1 / (val * L)
}
S7::method(d2linkfun, ClogLogLink) <- function(x, theta) {
  val <- 1 - theta
  L <- log(val)
  -(L + 1) / ((val^2) * (L^2))
}
S7::method(d3linkfun, ClogLogLink) <- function(x, theta) {
  val <- 1 - theta
  L <- log(val)
  -(2 * (L^2) + 3 * L + 2) / ((val^3) * (L^3))
}
S7::method(d4linkfun, ClogLogLink) <- function(x, theta) {
  val <- 1 - theta
  L <- log(val)
  -(6 * (L^3) + 11 * (L^2) + 12 * L + 6) / ((val^4) * (L^4))
}

# Exact analytical derivatives of the inverse link function (wrt eta).
#
# Written as T_k(eta) = exp(k*eta - exp(eta)), which sidesteps the "Inf * 0" NaN
# that comes of forming exp(k*eta) and exp(-exp(eta)) separately. The floor that
# used to be applied to z here was a no-op -- subtracting 2.2e-16 from eta cannot
# change exp(eta) -- so it is gone.
S7::method(dlinkinv, ClogLogLink) <- function(x, eta) {
  z <- exp(eta)
  exp(eta - z)
}
S7::method(d2linkinv, ClogLogLink) <- function(x, eta) {
  z <- exp(eta)
  exp(eta - z) - exp(2 * eta - z)
}
S7::method(d3linkinv, ClogLogLink) <- function(x, eta) {
  z <- exp(eta)
  exp(eta - z) - 3 * exp(2 * eta - z) + exp(3 * eta - z)
}
S7::method(d4linkinv, ClogLogLink) <- function(x, eta) {
  z <- exp(eta)
  exp(eta - z) - 7 * exp(2 * eta - z) + 6 * exp(3 * eta - z) - exp(4 * eta - z)
}

#' @title The Complementary Log-Log (ClogLog) Link Function
#'
#' @include generics.R
#' @include link_class.R
#' @description
#' Creates an S7 object of class \code{link} implementing the Complementary Log-Log 
#' transformation. This link is highly asymmetric and is predominantly used for modeling 
#' binary data where the probability of the event approaches 1 very slowly but approaches 
#' 0 rather sharply.
#'
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

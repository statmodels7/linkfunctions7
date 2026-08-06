#' @title S7 Class for the Softplus Link
#'
#' @description
#' The class \code{\link{softplus_link}} instantiates. Carries the scale
#' parameter \code{a} in a dedicated property.
#'
#' @param a The scale parameter, strictly positive.
#'
#' @return An S7 object of class \code{SoftplusLink}, inheriting from
#'   \code{\link{link}}.
#'
#' @seealso \code{\link{softplus_link}}, the constructor users call.
#' @keywords internal
SoftplusLink <- S7::new_class(
  name = "SoftplusLink",
  parent = link,
  properties = list(
    a = S7::class_numeric
  ),
  validator = function(self) {
    if (length(self@a) != 1L || is.na(self@a) || self@a <= 0) {
      "Scale parameter 'a' must be a single number strictly greater than 0."
    }
  }
)

# --- Methods for SoftplusLink ---
#
# Everything on the forward side is written in u = -expm1(-z), z = a * theta,
# rather than in val = expm1(z). The two are algebraically identical, but val
# overflows: expm1(z) is Inf for z > 709, so linkfun() returned Inf and every
# derivative NaN. Because the derivatives divide by val^k, the higher orders went
# first -- d4linkfun was already NaN at theta = 177 with a = 1, and at theta = 24
# with a = 30, which is an entirely ordinary parameter value for a steep softplus.
#
# The identity that removes the overflow is
#
#     log(expm1(z)) = z + log(-expm1(-z)),
#
# exact at both ends: as z -> 0 the right side tends to log(z), and as z grows
# log(-expm1(-z)) -> 0 so the whole thing tends to z, which is what the softplus
# is asymptotically. Writing e = exp(-z), the derivatives follow by dividing
# numerator and denominator of each old expression by the appropriate power of
# exp(z), which is what turns every growing quantity into a decaying one.

S7::method(linkfun, SoftplusLink) <- function(x, theta) {
  z <- x@a * theta
  (z + log(-expm1(-z))) / x@a
}
S7::method(linkinv, SoftplusLink) <- function(x, eta) {
  # "Log-Sum-Exp" algebraic trick to avoid computing exp() of large numbers
  pmax(0, eta) + log1p(exp(-abs(x@a * eta))) / x@a
}

# Exact analytical derivatives of the link function (wrt theta)
S7::method(dlinkfun, SoftplusLink) <- function(x, theta) lk_softplus_fwd_cpp(theta, x@a, 1L)
S7::method(d2linkfun, SoftplusLink) <- function(x, theta) lk_softplus_fwd_cpp(theta, x@a, 2L)
S7::method(d3linkfun, SoftplusLink) <- function(x, theta) lk_softplus_fwd_cpp(theta, x@a, 3L)
S7::method(d4linkfun, SoftplusLink) <- function(x, theta) lk_softplus_fwd_cpp(theta, x@a, 4L)

# Exact analytical derivatives of the inverse link function (wrt eta).
#
# The softplus is an antiderivative of the logistic, so its inverse-link
# derivatives are the logistic ones shifted a place: h' = sigma, h^(k+1) =
# a^k sigma^(k). They are the same four polynomials the logit and the doubly
# bounded link use, hence the shared helper.
S7::method(dlinkinv, SoftplusLink) <- function(x, eta) {
  stats::plogis(x@a * eta)
}
S7::method(d2linkinv, SoftplusLink) <- function(x, eta) {
  x@a * lk_logit_inv_cpp(x@a * eta, 1L)
}
S7::method(d3linkinv, SoftplusLink) <- function(x, eta) {
  (x@a^2) * lk_logit_inv_cpp(x@a * eta, 2L)
}
S7::method(d4linkinv, SoftplusLink) <- function(x, eta) {
  (x@a^3) * lk_logit_inv_cpp(x@a * eta, 3L)
}

#' @title The Softplus Link Function
#'
#' @include generics.R
#' @include link_class.R
#' @description
#' The softplus link with scale \eqn{a}: the inverse is
#' \eqn{\theta = \log(1 + e^{a\eta})/a}, a smooth approximation of
#' \eqn{\max(0, \eta)} that sharpens as \eqn{a} grows.
#' @param a A numeric value specifying the scaling parameter (smoothness/steepness). 
#' Must be strictly positive. Defaults to 1.
#'
#' @details
#' The Softplus link describes the relationship where the response parameter \eqn{\theta} 
#' is the Softplus of the linear predictor \eqn{\eta}.
#'
#' Mathematically:
#' \itemize{
#'   \item Inverse Link (Softplus): \eqn{\theta = \frac{1}{a} \log(1 + \exp(a \eta))}
#'   \item Link Function: \eqn{\eta = \frac{1}{a} \log(\exp(a \theta) - 1)}
#' }
#'
#' \strong{Behavior:}
#' For large negative \eqn{\eta}, \eqn{\theta \approx 0}.
#' For large positive \eqn{\eta}, \eqn{\theta \approx \eta} (linear behavior), whereas 
#' a Log link would imply \eqn{\theta = \exp(\eta)} (exponential behavior).
#'
#' \strong{Numerical Stability:}
#' Both directions are written so that no intermediate quantity grows with
#' \eqn{a\theta} or \eqn{a\eta}. The inverse link uses the log-sum-exp form, and
#' the forward link and its derivatives are expressed in
#' \eqn{u = 1 - e^{-a\theta}} rather than in \eqn{e^{a\theta} - 1}, which
#' overflows once \eqn{a\theta} passes about 709 — and, because the derivatives
#' divide by its fourth power, well before that at the higher orders.
#'
#' The mathematical domain of \eqn{\theta} is \code{c(0, Inf)}.
#'
#' @return An S7 object of class \code{SoftplusLink} (inheriting from \code{link}) containing the transformation functions,
#' their exact analytical derivatives up to the fourth order, and the parameter \code{a}.
#'
#' @examples
#' lk <- softplus_link(a = 2)
#' lk
#'
#' theta <- c(0.5, 1, 5)
#' eta <- linkfun(lk, theta)
#' eta
#' linkinv(lk, eta)          # back to theta
#'
#' # derivatives of either direction, to fourth order
#' dlinkfun(lk, theta)
#' d4linkinv(lk, eta)
#'
#' # unlike the log link, softplus is asymptotically linear in eta
#' linkinv(softplus_link(), c(1, 10, 100))
#'
#' @seealso \code{\link{link}}, \code{\link{log_link}}, \code{\link{identity_link}}
#' @importFrom stats plogis
#' @export
softplus_link <- function(a = 1) {
  # The class validator enforces this too; checking here as well is what makes
  # the message name the constructor the user actually called.
  if (length(a) != 1L || !is.numeric(a) || is.na(a) || a <= 0) {
    stop("Scale parameter 'a' must be a single number strictly greater than 0.",
         call. = FALSE)
  }

  SoftplusLink(
    link_name = paste0("softplus(a=", round(a, 5), ")"),
    link_bounds = c(0, Inf),
    link_params = list(a = a),
    a = a
  )
}
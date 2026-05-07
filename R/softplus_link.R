SoftplusLink <- S7::new_class(
  name = "SoftplusLink",
  parent = link,
  properties = list(
    a = S7::class_numeric
  ),
  validator = function(self) {
    if (self@a <= 0) {
      "Scale parameter 'a' must be strictly greater than 0."
    }
  }
)

# --- Methods for SoftplusLink ---

# Forward and inverse link functions with numerical stability guards
S7::method(linkfun, SoftplusLink) <- function(x, theta) {
  log(expm1(x@a * theta)) / x@a
}
S7::method(linkinv, SoftplusLink) <- function(x, eta) {
  # "Log-Sum-Exp" algebraic trick to avoid computing exp() of large numbers
  pmax(0, eta) + log(1 + exp(-abs(x@a * eta))) / x@a
}

# Exact analytical derivatives of the link function (wrt theta)
# We use val = expm1(a * theta) to express derivatives safely and efficiently
S7::method(dlinkfun, SoftplusLink) <- function(x, theta) {
  val <- expm1(x@a * theta); (val + 1) / val
}
S7::method(d2linkfun, SoftplusLink) <- function(x, theta) {
  val <- expm1(x@a * theta); -x@a * (val + 1) / (val^2)
}
S7::method(d3linkfun, SoftplusLink) <- function(x, theta) {
  val <- expm1(x@a * theta); (x@a^2) * (val + 1) * (val + 2) / (val^3)
}
S7::method(d4linkfun, SoftplusLink) <- function(x, theta) {
  val <- expm1(x@a * theta); -(x@a^3) * (val + 1) * (val^2 + 6 * val + 6) / (val^4)
}

# Exact analytical derivatives of the inverse link function (wrt eta)
# They cleanly resolve into polynomials of the logistic function p
S7::method(dlinkinv, SoftplusLink) <- function(x, eta) {
  stats::plogis(x@a * eta)
}
S7::method(d2linkinv, SoftplusLink) <- function(x, eta) {
  p <- stats::plogis(x@a * eta); x@a * p * (1 - p)
}
S7::method(d3linkinv, SoftplusLink) <- function(x, eta) {
  p <- stats::plogis(x@a * eta); (x@a^2) * p * (1 - p) * (1 - 2 * p)
}
S7::method(d4linkinv, SoftplusLink) <- function(x, eta) {
  p <- stats::plogis(x@a * eta); (x@a^3) * p * (1 - p) * (1 - 6 * p + 6 * (p^2))
}

#' @title The Softplus Link Function
#'
#' @include generics.R
#' @include link_class.R
#' @description
#' Creates an S7 object of class \code{link} implementing the Softplus transformation.
#' The Softplus function is a smooth approximation of the rectifier (ReLU) and ensures 
#' the parameter \eqn{\theta} remains positive. Unlike the Log link, which implies an 
#' exponential relationship throughout, the Softplus link asymptotically approaches 
#' linearity for large positive values of the linear predictor.
#'
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
#' The inverse link implementation intelligently utilizes conditional algebraic logic 
#' to ensure robust numerical stability for large positive values of \eqn{\eta}, 
#' entirely avoiding precision overflow.
#'
#' The mathematical domain of \eqn{\theta} is \code{c(0, Inf)}.
#'
#' @return An S7 object of class \code{SoftplusLink} (inheriting from \code{link}) containing the transformation functions,
#' their exact analytical derivatives up to the fourth order, and the parameter \code{a}.
#'
#' @seealso \code{\link{link}}, \code{\link{log_link}}, \code{\link{identity_link}}
#' @importFrom stats plogis
#' @export
softplus_link <- function(a = 1) {
  if (a <= 0) {
    stop("Scale parameter 'a' must be strictly greater than 0.")
  }

  SoftplusLink(
    link_name = paste0("softplus(a=", round(a, 5), ")"),
    link_bounds = c(0, Inf),
    link_params = list(a = a),
    a = a
  )
}
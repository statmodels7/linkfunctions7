ClogLogLink <- S7::new_class(
  name = "ClogLogLink",
  parent = link
)

# --- Methods for ClogLogLink ---

S7::method(linkfun, ClogLogLink) <- function(x, theta) {
  log(-log(1 - theta))
}

S7::method(linkinv, ClogLogLink) <- function(x, eta) {
  # Bounded strictly between Machine EPS and 1 - EPS to prevent catastrophic 0 or 1
  pmax(pmin(1 - exp(-exp(eta)), 1 - .Machine$double.eps), .Machine$double.eps)
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

# Exact analytical derivatives of the inverse link function (wrt eta)
# We utilize safe exponential forms T_k(eta) = exp(k*eta - exp(eta)). 
# This elegantly sidesteps NaN evaluations caused by "Inf * 0" during overflows.
S7::method(dlinkinv, ClogLogLink) <- function(x, eta) {
  z <- pmax(exp(eta), .Machine$double.eps)
  exp(eta - z)
}
S7::method(d2linkinv, ClogLogLink) <- function(x, eta) {
  z <- pmax(exp(eta), .Machine$double.eps)
  exp(eta - z) - exp(2 * eta - z)
}
S7::method(d3linkinv, ClogLogLink) <- function(x, eta) {
  z <- pmax(exp(eta), .Machine$double.eps)
  exp(eta - z) - 3 * exp(2 * eta - z) + exp(3 * eta - z)
}
S7::method(d4linkinv, ClogLogLink) <- function(x, eta) {
  z <- pmax(exp(eta), .Machine$double.eps)
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
#' and their exact analytical derivatives up to the fourth order.
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
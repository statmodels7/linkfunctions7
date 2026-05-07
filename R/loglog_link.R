LogLogLink <- S7::new_class(
  name = "LogLogLink",
  parent = link
)

# --- Methods for LogLogLink ---

S7::method(linkfun, LogLogLink) <- function(x, theta) -log(-log(theta))
S7::method(linkinv, LogLogLink) <- function(x, eta) exp(-exp(-eta))

# Exact analytical derivatives of the link function (wrt theta)
S7::method(dlinkfun, LogLogLink) <- function(x, theta) {
  -1 / (theta * log(theta))
}
S7::method(d2linkfun, LogLogLink) <- function(x, theta) {
  l <- log(theta)
  (1 + l) / (theta^2 * (l^2))
}
S7::method(d3linkfun, LogLogLink) <- function(x, theta) {
  l <- log(theta)
  -(2 + 3 * l + 2 * (l^2)) / (theta^3 * (l^3))
}
S7::method(d4linkfun, LogLogLink) <- function(x, theta) {
  l <- log(theta)
  (6 + 12 * l + 11 * (l^2) + 6 * (l^3)) / (theta^4 * (l^4))
}

# Exact analytical derivatives of the inverse link function (wrt eta)
# Utilizing the term z = exp(-eta) to evaluate derivatives as polynomials,
# thus maximizing computational performance.
S7::method(dlinkinv, LogLogLink) <- function(x, eta) { z <- exp(-eta); exp(-z) * z }
S7::method(d2linkinv, LogLogLink) <- function(x, eta) { z <- exp(-eta); exp(-z) * (z^2 - z) }
S7::method(d3linkinv, LogLogLink) <- function(x, eta) { z <- exp(-eta); exp(-z) * (z^3 - 3 * z^2 + z) }
S7::method(d4linkinv, LogLogLink) <- function(x, eta) { z <- exp(-eta); exp(-z) * (z^4 - 6 * z^3 + 7 * z^2 - z) }

#' @title The Log-Log Link Function
#'
#' @include generics.R
#' @include link_class.R
#' @description
#' Creates an S7 object of class \code{link} implementing the Log-Log transformation.
#' This link is the asymmetric opposite of the \code{cloglog} link and is often used 
#' in survival analysis (e.g., Gompertz models) or for modeling binary data with 
#' specific asymmetry requirements.
#'
#' @details
#' The Log-Log link is mathematically defined as \eqn{\eta = -\log(-\log(\theta))}.
#' Consequently, the inverse link is derived as \eqn{\theta = \exp(-\exp(-\eta))}.
#'
#' \strong{Asymmetry:} Unlike the symmetric Logit or Probit links, the Log-Log link 
#' is highly asymmetric. It is particularly suitable for modeling events where the 
#' probability approaches 0 very slowly but approaches 1 sharply. This represents 
#' the mathematical reverse of the \code{cloglog} link.
#'
#' The strict mathematical domain for \eqn{\theta} is \code{c(0, 1)}.
#'
#' @return An S7 object of class \code{LogLogLink} (inheriting from \code{link}) containing the transformation functions
#' and their exact analytical derivatives up to the fourth order.
#'
#' @seealso \code{\link{link}}, \code{\link{cloglog_link}}, \code{\link{logit_link}}
#' @export
loglog_link <- function() {
  LogLogLink(
    link_name = "loglog",
    link_bounds = c(0, 1),
    link_params = NULL
  )
}
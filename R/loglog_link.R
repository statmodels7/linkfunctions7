#' @title S7 Class for the LogLog Link
#'
#' @description
#' The class \code{\link{loglog_link}} instantiates.
#'
#' @return An S7 object of class \code{LogLogLink}, inheriting from
#'   \code{\link{link}}.
#'
#' @seealso \code{\link{loglog_link}}, the constructor users call.
#' @keywords internal
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
#' The log-log link \eqn{\eta = -\log(-\log\theta)} on \eqn{(0, 1)},
#' with inverse \eqn{\theta = \exp(-e^{-\eta})}; the mirror image of
#' \code{\link{cloglog_link}}.
#' @details
#' The Log-Log link is mathematically defined as \eqn{\eta = -\log(-\log(\theta))}.
#' Consequently, the inverse link is derived as \eqn{\theta = \exp(-\exp(-\eta))}.
#'
#' Unlike the logit and the probit the link is asymmetric: the probability
#' approaches 0 slowly and 1 sharply, the mirror image of
#' \code{\link{cloglog_link}}. The domain of \eqn{\theta} is \eqn{(0, 1)}.
#'
#' @return An S7 object of class \code{LogLogLink} (inheriting from \code{link}) containing the transformation functions
#' and their exact analytical derivatives up to the fourth order.
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
#' @seealso \code{\link{link}}, \code{\link{cloglog_link}}, \code{\link{logit_link}}
#' @export
loglog_link <- function() {
  LogLogLink(
    link_name = "loglog",
    link_bounds = c(0, 1),
    link_params = NULL
  )
}

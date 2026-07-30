#' @title S7 Class for the Logit Link
#'
#' @description
#' The class \code{\link{logit_link}} instantiates.
#'
#' @return An S7 object of class \code{LogitLink}, inheriting from
#'   \code{\link{link}}.
#'
#' @seealso \code{\link{logit_link}}, the constructor users call.
#' @keywords internal
LogitLink <- S7::new_class(
  name = "LogitLink",
  parent = link
)

# --- Methods for LogitLink ---

# Forward and inverse link functions using native R C-level functions
S7::method(linkfun, LogitLink) <- function(x, theta) stats::qlogis(theta)
S7::method(linkinv, LogitLink) <- function(x, eta) stats::plogis(eta)

# Exact analytical derivatives of the link function (wrt theta)
S7::method(dlinkfun, LogitLink) <- function(x, theta) 1 / (theta * (1 - theta))
S7::method(d2linkfun, LogitLink) <- function(x, theta) (2 * theta - 1) / ((theta * (1 - theta))^2)
S7::method(d3linkfun, LogitLink) <- function(x, theta) 2 / (theta^3) + 2 / ((1 - theta)^3)
S7::method(d4linkfun, LogitLink) <- function(x, theta) -6 / (theta^4) + 6 / ((1 - theta)^4)

# Exact analytical derivatives of the inverse link function (wrt eta).
#
# Polynomials in the probability p itself. They are shared with the doubly
# bounded link, which scales them by the interval width, and with the softplus,
# which uses them one order down; see logistic_deriv().
S7::method(dlinkinv, LogitLink) <- function(x, eta) {
  logistic_deriv(stats::plogis(eta), 1L)
}
S7::method(d2linkinv, LogitLink) <- function(x, eta) {
  logistic_deriv(stats::plogis(eta), 2L)
}
S7::method(d3linkinv, LogitLink) <- function(x, eta) {
  logistic_deriv(stats::plogis(eta), 3L)
}
S7::method(d4linkinv, LogitLink) <- function(x, eta) {
  logistic_deriv(stats::plogis(eta), 4L)
}

#' @title The Logit Link Function
#'
#' @include generics.R
#' @include link_class.R
#' @description
#' Creates an S7 object of class \code{link} implementing the Logit (log-odds) transformation.
#' This is the canonical link function for the success probability parameter of the Bernoulli 
#' and Binomial distributions and serves as the foundation of Logistic Regression.
#'
#' @details
#' The Logit link is defined mathematically as \eqn{\eta = \log(\frac{\theta}{1 - \theta})}.
#' The inverse link is the standard logistic function (sigmoid): 
#' \eqn{\theta = \frac{1}{1 + \exp(-\eta)}}.
#'
#' The link is perfectly symmetric around \eqn{\theta = 0.5} (which corresponds to \eqn{\eta = 0}). 
#' It elegantly interprets the linear predictor as the log-odds of the event probability.
#'
#' The mathematical domain of \eqn{\theta} is \code{c(0, 1)}.
#'
#' \strong{Implementation Details:} 
#' This function internally delegates to R's highly optimized native functions 
#' \code{stats::qlogis} and \code{stats::plogis}. This ensures maximum numerical stability 
#' and precision, especially when evaluating probabilities exceedingly close to the 
#' boundaries of 0 and 1.
#'
#' @return An S7 object of class \code{LogitLink} (inheriting from \code{link}) containing the transformation functions
#' and their exact analytical derivatives up to the fourth order.
#'
#' @examples
#' lk <- logit_link()
#' lk
#'
#' p <- c(0.1, 0.5, 0.9)
#' eta <- linkfun(lk, p)      # log-odds
#' eta
#' linkinv(lk, eta)
#'
#' # the inverse derivatives are polynomials in p: the first is the variance
#' # of a Bernoulli, p(1 - p)
#' dlinkinv(lk, 0)
#'
#' # all four orders at once
#' vapply(1:4, function(k) linkinvderiv(lk, 0, order = k), numeric(1))
#'
#' # every mathematical property is checkable
#' check_link(lk)
#'
#' @seealso \code{\link{link}}, \code{\link{probit_link}}, \code{\link{cloglog_link}}
#' @importFrom stats qlogis plogis
#' @export
logit_link <- function() {
  LogitLink(
    link_name = "logit",
    link_bounds = c(0, 1),
    
    # The logit link requires no additional mathematical parameters
    link_params = NULL
  )
}
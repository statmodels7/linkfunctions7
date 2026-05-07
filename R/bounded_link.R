DoublyBoundedLink <- S7::new_class(
  name = "DoublyBoundedLink",
  parent = link,
  properties = list(
    lwr = S7::class_numeric,
    upr = S7::class_numeric
  ),
  validator = function(self) {
    if (self@lwr >= self@upr) {
      "Lower bound 'lwr' must be strictly less than upper bound 'upr'."
    }
  }
)

LowerBoundedLink <- S7::new_class(
  name = "LowerBoundedLink",
  parent = link,
  properties = list(
    lwr = S7::class_numeric
  )
)

UpperBoundedLink <- S7::new_class(
  name = "UpperBoundedLink",
  parent = link,
  properties = list(
    upr = S7::class_numeric
  )
)

# --- Methods for DoublyBoundedLink ---

S7::method(linkfun, DoublyBoundedLink) <- function(x, theta) {
  W <- x@upr - x@lwr
  p <- (theta - x@lwr) / W
  stats::qlogis(p)
}
S7::method(linkinv, DoublyBoundedLink) <- function(x, eta) {
  W <- x@upr - x@lwr
  x@lwr + W * stats::plogis(eta)
}
S7::method(dlinkfun, DoublyBoundedLink) <- function(x, theta) {
  W <- x@upr - x@lwr
  p <- (theta - x@lwr) / W
  1 / (W * p * (1 - p))
}
S7::method(d2linkfun, DoublyBoundedLink) <- function(x, theta) {
  W <- x@upr - x@lwr
  p <- (theta - x@lwr) / W
  (2 * p - 1) / ((W^2) * (p^2) * ((1 - p)^2))
}
S7::method(d3linkfun, DoublyBoundedLink) <- function(x, theta) {
  W <- x@upr - x@lwr
  p <- (theta - x@lwr) / W
  (2 / (p^3) + 2 / ((1 - p)^3)) / (W^3)
}
S7::method(d4linkfun, DoublyBoundedLink) <- function(x, theta) {
  W <- x@upr - x@lwr
  p <- (theta - x@lwr) / W
  (-6 / (p^4) + 6 / ((1 - p)^4)) / (W^4)
}
S7::method(dlinkinv, DoublyBoundedLink) <- function(x, eta) {
  p <- stats::plogis(eta)
  (x@upr - x@lwr) * p * (1 - p)
}
S7::method(d2linkinv, DoublyBoundedLink) <- function(x, eta) {
  p <- stats::plogis(eta)
  (x@upr - x@lwr) * p * (1 - p) * (1 - 2 * p)
}
S7::method(d3linkinv, DoublyBoundedLink) <- function(x, eta) {
  p <- stats::plogis(eta)
  (x@upr - x@lwr) * p * (1 - p) * (1 - 6 * p + 6 * (p^2))
}
S7::method(d4linkinv, DoublyBoundedLink) <- function(x, eta) {
  p <- stats::plogis(eta)
  (x@upr - x@lwr) * p * (1 - p) * (1 - 14 * p + 36 * (p^2) - 24 * (p^3))
}

# --- Methods for LowerBoundedLink ---

S7::method(linkfun, LowerBoundedLink) <- function(x, theta) log(theta - x@lwr)
S7::method(linkinv, LowerBoundedLink) <- function(x, eta) x@lwr + pmax(exp(eta), .Machine$double.eps)
S7::method(dlinkfun, LowerBoundedLink) <- function(x, theta) 1 / (theta - x@lwr)
S7::method(d2linkfun, LowerBoundedLink) <- function(x, theta) -1 / ((theta - x@lwr)^2)
S7::method(d3linkfun, LowerBoundedLink) <- function(x, theta) 2 / ((theta - x@lwr)^3)
S7::method(d4linkfun, LowerBoundedLink) <- function(x, theta) -6 / ((theta - x@lwr)^4)
S7::method(dlinkinv, LowerBoundedLink) <- function(x, eta) pmax(exp(eta), .Machine$double.eps)
S7::method(d2linkinv, LowerBoundedLink) <- function(x, eta) pmax(exp(eta), .Machine$double.eps)
S7::method(d3linkinv, LowerBoundedLink) <- function(x, eta) pmax(exp(eta), .Machine$double.eps)
S7::method(d4linkinv, LowerBoundedLink) <- function(x, eta) pmax(exp(eta), .Machine$double.eps)

# --- Methods for UpperBoundedLink ---

S7::method(linkfun, UpperBoundedLink) <- function(x, theta) log(x@upr - theta)
S7::method(linkinv, UpperBoundedLink) <- function(x, eta) x@upr - pmax(exp(eta), .Machine$double.eps)
S7::method(dlinkfun, UpperBoundedLink) <- function(x, theta) -1 / (x@upr - theta)
S7::method(d2linkfun, UpperBoundedLink) <- function(x, theta) -1 / ((x@upr - theta)^2)
S7::method(d3linkfun, UpperBoundedLink) <- function(x, theta) -2 / ((x@upr - theta)^3)
S7::method(d4linkfun, UpperBoundedLink) <- function(x, theta) -6 / ((x@upr - theta)^4)
S7::method(dlinkinv, UpperBoundedLink) <- function(x, eta) -pmax(exp(eta), .Machine$double.eps)
S7::method(d2linkinv, UpperBoundedLink) <- function(x, eta) -pmax(exp(eta), .Machine$double.eps)
S7::method(d3linkinv, UpperBoundedLink) <- function(x, eta) -pmax(exp(eta), .Machine$double.eps)
S7::method(d4linkinv, UpperBoundedLink) <- function(x, eta) -pmax(exp(eta), .Machine$double.eps)

# --- The General Bounded Link Function Factory ---

#' @title The General Bounded Link Function
#'
#' @include generics.R
#' @include link_class.R
#' @description
#' Creates an S7 object of class \code{link} that maps a constrained interval to the 
#' whole real line. By specifying \code{lwr} and \code{upr}, this function dynamically 
#' constructs a doubly bounded (scaled logit), lower bounded (shifted log), upper bounded, 
#' or unbounded (identity) link function.
#'
#' @param lwr Numeric or \code{NULL}. The lower bound of the interval.
#' @param upr Numeric or \code{NULL}. The upper bound of the interval.
#'
#' @details
#' \strong{Doubly Bounded (\code{lwr} and \code{upr} provided):}
#' Transforms \eqn{\theta} by normalizing it to \code{c(0, 1)} via 
#' \eqn{p = \frac{\theta - \text{lwr}}{\text{upr} - \text{lwr}}}, and then applying the logit function.
#' 
#' \strong{Lower Bounded (\code{lwr} provided, \code{upr = NULL}):}
#' Defined as \eqn{\eta = \log(\theta - \text{lwr})}, with inverse \eqn{\theta = \exp(\eta) + \text{lwr}}.
#' 
#' \strong{Upper Bounded (\code{lwr = NULL}, \code{upr} provided):}
#' Defined as \eqn{\eta = \log(\text{upr} - \theta)}, with inverse \eqn{\theta = \text{upr} - \exp(\eta)}.
#' 
#' \strong{Unbounded (\code{lwr = NULL}, \code{upr = NULL}):}
#' Returns the standard \code{\link{identity_link}}.
#'
#' @return An S7 object of class \code{link} containing the transformation functions
#' and their exact analytical derivatives up to the fourth order.
#'
#' @seealso \code{\link{link}}, \code{\link{logit_link}}, \code{\link{log_link}}
#' @importFrom stats qlogis plogis
#' @export
bounded_link <- function(lwr = NULL, upr = NULL) {
  
  # Case 0: Unbounded link effectively reduces to Identity
  if (is.null(lwr) && is.null(upr)) {
    return(identity_link())
  }
  
  # Case 1: Doubly Bounded
  if (!is.null(lwr) && !is.null(upr)) {
    return(DoublyBoundedLink(
      link_name = paste0("bounded(lwr=", lwr, ", upr=", upr, ")"),
      link_bounds = c(lwr, upr),
      link_params = list(lwr = lwr, upr = upr),
      lwr = lwr,
      upr = upr
    ))
  }
  
  # Case 2: Lower Bounded
  if (!is.null(lwr) && is.null(upr)) {
    return(LowerBoundedLink(
      link_name = paste0("lower_bounded(lwr=", lwr, ")"),
      link_bounds = c(lwr, Inf),
      link_params = list(lwr = lwr),
      lwr = lwr
    ))
  }
  
  # Case 3: Upper Bounded
  if (is.null(lwr) && !is.null(upr)) {
    return(UpperBoundedLink(
      link_name = paste0("upper_bounded(upr=", upr, ")"),
      link_bounds = c(-Inf, upr),
      link_params = list(upr = upr),
      upr = upr
    ))
  }
}
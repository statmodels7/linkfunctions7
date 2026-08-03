#' @title S7 Class for a Doubly Bounded Link
#'
#' @description
#' The class \code{\link{bounded_link}} instantiates when given both endpoints.
#' It is the logit of \eqn{p = (\theta - \mathrm{lwr})/W}, the position of
#' \eqn{\theta} within the interval.
#'
#' @details
#' The interval width \eqn{W = \mathrm{upr} - \mathrm{lwr}} is stored as its own
#' property rather than recomputed. Every method needs it, and reading two S7
#' properties and subtracting cost about a third of a call to \code{dlinkinv()};
#' the constructor is the only place it can change.
#'
#' @param lwr,upr The interval endpoints.
#' @param width The interval width, \code{upr - lwr}, set by the constructor.
#'
#' @return An S7 object of class \code{DoublyBoundedLink}, inheriting from
#'   \code{\link{link}}.
#'
#' @seealso \code{\link{bounded_link}}, the constructor users call.
#' @keywords internal
DoublyBoundedLink <- S7::new_class(
  name = "DoublyBoundedLink",
  parent = link,
  properties = list(
    lwr = S7::class_numeric,
    upr = S7::class_numeric,
    width = S7::class_numeric
  ),
  validator = function(self) {
    if (self@lwr >= self@upr) {
      "Lower bound 'lwr' must be strictly less than upper bound 'upr'."
    }
  }
)

#' @title S7 Class for a Lower Bounded Link
#'
#' @description
#' The class \code{\link{bounded_link}} instantiates when given only a lower
#' endpoint. It is the log link shifted to start at \code{lwr}.
#'
#' @param lwr The lower endpoint.
#'
#' @return An S7 object of class \code{LowerBoundedLink}, inheriting from
#'   \code{\link{link}}.
#'
#' @seealso \code{\link{bounded_link}}, the constructor users call.
#' @keywords internal
LowerBoundedLink <- S7::new_class(
  name = "LowerBoundedLink",
  parent = link,
  properties = list(
    lwr = S7::class_numeric
  )
)

#' @title S7 Class for an Upper Bounded Link
#'
#' @description
#' The class \code{\link{bounded_link}} instantiates when given only an upper
#' endpoint. It is the mirror image of \code{\link{LowerBoundedLink}}, the log of
#' the distance below \code{upr}.
#'
#' @param upr The upper endpoint.
#'
#' @return An S7 object of class \code{UpperBoundedLink}, inheriting from
#'   \code{\link{link}}.
#'
#' @seealso \code{\link{bounded_link}}, the constructor users call.
#' @keywords internal
UpperBoundedLink <- S7::new_class(
  name = "UpperBoundedLink",
  parent = link,
  properties = list(
    upr = S7::class_numeric
  )
)

# --- Methods for DoublyBoundedLink ---
#
# The forward derivatives are the logit's, divided by W^k because p is theta
# rescaled by W; the inverse ones are the logit's multiplied by W. Both sets are
# therefore the shared logistic polynomials rather than four more transcriptions.

S7::method(linkfun, DoublyBoundedLink) <- function(x, theta) {
  stats::qlogis((theta - x@lwr) / x@width)
}
S7::method(linkinv, DoublyBoundedLink) <- function(x, eta) {
  x@lwr + x@width * stats::plogis(eta)
}
S7::method(dlinkfun, DoublyBoundedLink) <- function(x, theta) {
  W <- x@width
  p <- (theta - x@lwr) / W
  1 / (W * p * (1 - p))
}
S7::method(d2linkfun, DoublyBoundedLink) <- function(x, theta) {
  W <- x@width
  p <- (theta - x@lwr) / W
  (2 * p - 1) / ((W^2) * (p^2) * ((1 - p)^2))
}
S7::method(d3linkfun, DoublyBoundedLink) <- function(x, theta) {
  W <- x@width
  p <- (theta - x@lwr) / W
  (2 / (p^3) + 2 / ((1 - p)^3)) / (W^3)
}
S7::method(d4linkfun, DoublyBoundedLink) <- function(x, theta) {
  W <- x@width
  p <- (theta - x@lwr) / W
  (-6 / (p^4) + 6 / ((1 - p)^4)) / (W^4)
}
S7::method(dlinkinv, DoublyBoundedLink) <- function(x, eta) {
  x@width * logistic_deriv(stats::plogis(eta), 1L)
}
S7::method(d2linkinv, DoublyBoundedLink) <- function(x, eta) {
  x@width * logistic_deriv(stats::plogis(eta), 2L)
}
S7::method(d3linkinv, DoublyBoundedLink) <- function(x, eta) {
  x@width * logistic_deriv(stats::plogis(eta), 3L)
}
S7::method(d4linkinv, DoublyBoundedLink) <- function(x, eta) {
  x@width * logistic_deriv(stats::plogis(eta), 4L)
}

# --- Methods for LowerBoundedLink ---
#
# The log link on theta - lwr. Every inverse derivative is exp(eta), floored for
# the same reason as in LogLink.

S7::method(linkfun, LowerBoundedLink) <- function(x, theta) log(theta - x@lwr)
S7::method(linkinv, LowerBoundedLink) <- function(x, eta) x@lwr + exp_floored(eta)
S7::method(dlinkfun, LowerBoundedLink) <- function(x, theta) 1 / (theta - x@lwr)
S7::method(d2linkfun, LowerBoundedLink) <- function(x, theta) -1 / ((theta - x@lwr)^2)
S7::method(d3linkfun, LowerBoundedLink) <- function(x, theta) 2 / ((theta - x@lwr)^3)
S7::method(d4linkfun, LowerBoundedLink) <- function(x, theta) -6 / ((theta - x@lwr)^4)
S7::method(dlinkinv, LowerBoundedLink) <- function(x, eta) exp_floored(eta)
S7::method(d2linkinv, LowerBoundedLink) <- function(x, eta) exp_floored(eta)
S7::method(d3linkinv, LowerBoundedLink) <- function(x, eta) exp_floored(eta)
S7::method(d4linkinv, LowerBoundedLink) <- function(x, eta) exp_floored(eta)

# --- Methods for UpperBoundedLink ---
#
# The mirror image: theta = upr - exp(eta), so every inverse derivative picks up
# a sign and the forward ones are the log link's in upr - theta.

S7::method(linkfun, UpperBoundedLink) <- function(x, theta) log(x@upr - theta)
S7::method(linkinv, UpperBoundedLink) <- function(x, eta) x@upr - exp_floored(eta)
S7::method(dlinkfun, UpperBoundedLink) <- function(x, theta) -1 / (x@upr - theta)
S7::method(d2linkfun, UpperBoundedLink) <- function(x, theta) -1 / ((x@upr - theta)^2)
S7::method(d3linkfun, UpperBoundedLink) <- function(x, theta) -2 / ((x@upr - theta)^3)
S7::method(d4linkfun, UpperBoundedLink) <- function(x, theta) -6 / ((x@upr - theta)^4)
S7::method(dlinkinv, UpperBoundedLink) <- function(x, eta) -exp_floored(eta)
S7::method(d2linkinv, UpperBoundedLink) <- function(x, eta) -exp_floored(eta)
S7::method(d3linkinv, UpperBoundedLink) <- function(x, eta) -exp_floored(eta)
S7::method(d4linkinv, UpperBoundedLink) <- function(x, eta) -exp_floored(eta)

# --- The General Bounded Link Function Factory ---

#' @title The General Bounded Link Function
#'
#' @include generics.R
#' @include link_class.R
#' @description
#' The link for a parameter confined to \eqn{(lwr, upr)}: a scaled logit when
#' both endpoints are finite, a shifted log \eqn{\eta = \log(\theta - lwr)}
#' when only the lower is, its mirror image \eqn{\eta = \log(upr - \theta)}
#' when only the upper is, and the identity when neither is given.
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
#' and their exact analytical derivatives up to the fourth order. Which class
#' exactly depends on the endpoints given: \code{\link{DoublyBoundedLink}},
#' \code{\link{LowerBoundedLink}}, \code{\link{UpperBoundedLink}}, or an
#' \code{\link{IdentityLink}} when neither endpoint is supplied.
#'
#' @examples
#' # a parameter known to lie in (0, 10)
#' lk <- bounded_link(lwr = 0, upr = 10)
#' lk
#' linkinv(lk, c(-2, 0, 2))     # always inside the interval
#' linkfun(lk, c(1, 5, 9))
#'
#' # one-sided: a variance component bounded below by zero
#' bounded_link(lwr = 0)
#'
#' # no endpoints at all is the identity
#' bounded_link()
#'
#' # derivatives, as for any other link
#' dlinkinv(lk, 0)
#'
#' @seealso \code{\link{link}}, \code{\link{logit_link}}, \code{\link{log_link}}
#' @importFrom stats qlogis plogis
#' @export
bounded_link <- function(lwr = NULL, upr = NULL) {

  # A bound that is not a single finite number produces a link object whose
  # methods return NaN everywhere, which is a much harder thing to diagnose from
  # a model fit than a message here.
  ok <- function(v, nm) {
    if (length(v) != 1L || !is.numeric(v) || !is.finite(v)) {
      stop("'", nm, "' must be a single finite number.", call. = FALSE)
    }
  }
  if (!is.null(lwr)) ok(lwr, "lwr")
  if (!is.null(upr)) ok(upr, "upr")

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
      upr = upr,
      width = upr - lwr
    ))
  }

  # Case 2: Lower Bounded
  if (!is.null(lwr)) {
    return(LowerBoundedLink(
      link_name = paste0("lower_bounded(lwr=", lwr, ")"),
      link_bounds = c(lwr, Inf),
      link_params = list(lwr = lwr),
      lwr = lwr
    ))
  }

  # Case 3: Upper Bounded -- the only case left
  UpperBoundedLink(
    link_name = paste0("upper_bounded(upr=", upr, ")"),
    link_bounds = c(-Inf, upr),
    link_params = list(upr = upr),
    upr = upr
  )
}
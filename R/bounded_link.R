#' @title S7 Class for a Doubly Bounded Link
#'
#' @description
#' Carries the scaled logit on a finite interval \eqn{(l, u)}: the parameter is
#' mapped to its position \eqn{p = (\theta - l)/(u - l)} within the interval and
#' that position is carried by the logit, so
#' \eqn{\eta = \log(p/(1-p))}.
#'
#' Its derivatives are the logistic polynomials scaled by the width
#' \eqn{u - l}, so a bounded link costs no more than a logit.
#'
#' @details
#' The interval width \eqn{W = \mathrm{upr} - \mathrm{lwr}} is stored as its
#' own property, set by the constructor, so that each method reads one
#' property instead of two.
#'
#' @param lwr,upr The interval endpoints.
#' @param link_name A character string naming the link, set by the
#'   constructor and shown by `print()`.
#' @param link_bounds A length-two numeric vector, the open interval in
#'   which the parameter lies. Set by the constructor; see Value for this link's.
#' @param link_params A list of the link's own parameters, empty where it has
#'   none. Set by the constructor.
#' @param width The interval width, `upr - lwr`, set by the constructor.
#'
#' @return An S7 object of class `DoublyBoundedLink`, inheriting from [link()] and
#'   carrying its three properties `link_name`, `link_bounds` and
#'   `link_params`. Its `link_bounds` are the endpoints given, and its `link_params` holds
#'   `lwr` and `upr`.
#'
#' @section Methods:
#' Twelve methods are registered on this class: [linkfun()] and [linkinv()],
#' and the five derivative orders in each direction, [dlinkfun()] through
#' [d5linkfun()] going out and [dlinkinv()] through [d5linkinv()] coming
#' back. The forward derivatives are the logit's divided by \eqn{W^k} and
#' the inverse ones the logit's multiplied by \eqn{W}, because \eqn{p} is
#' \eqn{\theta} rescaled by the width. The inverse set therefore calls the
#' compiled kernel of the logit, and the forward set is the logit's
#' expressions written in \eqn{p}.
#'
#' @aliases linkfun.DoublyBoundedLink
#' @aliases linkinv.DoublyBoundedLink
#' @aliases dlinkfun.DoublyBoundedLink
#' @aliases d2linkfun.DoublyBoundedLink
#' @aliases d3linkfun.DoublyBoundedLink
#' @aliases d4linkfun.DoublyBoundedLink
#' @aliases d5linkfun.DoublyBoundedLink
#' @aliases dlinkinv.DoublyBoundedLink
#' @aliases d2linkinv.DoublyBoundedLink
#' @aliases d3linkinv.DoublyBoundedLink
#' @aliases d4linkinv.DoublyBoundedLink
#' @aliases d5linkinv.DoublyBoundedLink
#'
#' @seealso [bounded_link()], the constructor users call.
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
#' Carries the shifted log on \eqn{(l, \infty)}: \eqn{\eta = \log(\theta - l)},
#' with inverse \eqn{\theta = l + e^{\eta}}.
#'
#' Every derivative is the log link's, the shift being a constant that
#' differentiates away, and the exponential is floored at [exp_floor()] so the
#' parameter never reaches the bound exactly.
#'
#' @param link_name A character string naming the link, set by the
#'   constructor and shown by `print()`.
#' @param link_bounds A length-two numeric vector, the open interval in
#'   which the parameter lies. Set by the constructor; see Value for this link's.
#' @param link_params A list of the link's own parameters, empty where it has
#'   none. Set by the constructor.
#' @param lwr The lower endpoint.
#'
#' @return An S7 object of class `LowerBoundedLink`, inheriting from [link()] and
#'   carrying its three properties `link_name`, `link_bounds` and
#'   `link_params`. Its `link_bounds` are `c(lwr, Inf)`, and its `link_params` holds `lwr`.
#'
#' @section Methods:
#' Twelve methods are registered on this class: [linkfun()] and [linkinv()],
#' and the five derivative orders in each direction, [dlinkfun()] through
#' [d5linkfun()] going out and [dlinkinv()] through [d5linkinv()] coming
#' back. The forward derivatives are the log link's read at \eqn{\theta -
#' \mathrm{lwr}}, the shift being a constant that differentiates away. Every
#' inverse derivative is [exp_floored()] of \eqn{\eta} itself, the
#' exponential being its own derivative to every order.
#'
#' @aliases linkfun.LowerBoundedLink
#' @aliases linkinv.LowerBoundedLink
#' @aliases dlinkfun.LowerBoundedLink
#' @aliases d2linkfun.LowerBoundedLink
#' @aliases d3linkfun.LowerBoundedLink
#' @aliases d4linkfun.LowerBoundedLink
#' @aliases d5linkfun.LowerBoundedLink
#' @aliases dlinkinv.LowerBoundedLink
#' @aliases d2linkinv.LowerBoundedLink
#' @aliases d3linkinv.LowerBoundedLink
#' @aliases d4linkinv.LowerBoundedLink
#' @aliases d5linkinv.LowerBoundedLink
#'
#' @seealso [bounded_link()], the constructor users call.
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
#' Carries the reflected log on \eqn{(-\infty, u)}: \eqn{\eta = \log(u - \theta)},
#' with inverse \eqn{\theta = u - e^{\eta}}.
#'
#' It is the mirror image of [LowerBoundedLink()], the log of the distance
#' below `upr`. The reflection makes the map decreasing, so the odd-order
#' forward derivatives change sign relative to those of [LowerBoundedLink]
#' while the even ones do not, and every inverse derivative changes sign.
#'
#' @param link_name A character string naming the link, set by the
#'   constructor and shown by `print()`.
#' @param link_bounds A length-two numeric vector, the open interval in
#'   which the parameter lies. Set by the constructor; see Value for this link's.
#' @param link_params A list of the link's own parameters, empty where it has
#'   none. Set by the constructor.
#' @param upr The upper endpoint.
#'
#' @return An S7 object of class `UpperBoundedLink`, inheriting from [link()] and
#'   carrying its three properties `link_name`, `link_bounds` and
#'   `link_params`. Its `link_bounds` are `c(-Inf, upr)`, and its `link_params` holds `upr`.
#'
#' @section Methods:
#' Twelve methods are registered on this class: [linkfun()] and [linkinv()],
#' and the five derivative orders in each direction, [dlinkfun()] through
#' [d5linkfun()] going out and [dlinkinv()] through [d5linkinv()] coming
#' back. The reflection makes the map decreasing, so every inverse
#' derivative is the negative of [exp_floored()] of \eqn{\eta}, and the
#' forward ones
#' are the log's read at \eqn{\mathrm{upr} - \theta}, with the change of sign
#' that the reflection introduces at the odd orders.
#'
#' @aliases linkfun.UpperBoundedLink
#' @aliases linkinv.UpperBoundedLink
#' @aliases dlinkfun.UpperBoundedLink
#' @aliases d2linkfun.UpperBoundedLink
#' @aliases d3linkfun.UpperBoundedLink
#' @aliases d4linkfun.UpperBoundedLink
#' @aliases d5linkfun.UpperBoundedLink
#' @aliases dlinkinv.UpperBoundedLink
#' @aliases d2linkinv.UpperBoundedLink
#' @aliases d3linkinv.UpperBoundedLink
#' @aliases d4linkinv.UpperBoundedLink
#' @aliases d5linkinv.UpperBoundedLink
#'
#' @seealso [bounded_link()], the constructor users call.
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
S7::method(d5linkfun, DoublyBoundedLink) <- function(x, theta) {
  W <- x@width
  p <- (theta - x@lwr) / W
  (24 / (p^5) + 24 / ((1 - p)^5)) / (W^5)
}
S7::method(dlinkinv, DoublyBoundedLink) <- function(x, eta) {
  x@width * lk_logit_inv_cpp(eta, 1L)
}
S7::method(d2linkinv, DoublyBoundedLink) <- function(x, eta) {
  x@width * lk_logit_inv_cpp(eta, 2L)
}
S7::method(d3linkinv, DoublyBoundedLink) <- function(x, eta) {
  x@width * lk_logit_inv_cpp(eta, 3L)
}
S7::method(d4linkinv, DoublyBoundedLink) <- function(x, eta) {
  x@width * lk_logit_inv_cpp(eta, 4L)
}
S7::method(d5linkinv, DoublyBoundedLink) <- function(x, eta) {
  x@width * lk_logit_inv_cpp(eta, 5L)
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
S7::method(d5linkfun, LowerBoundedLink) <- function(x, theta) 24 / ((theta - x@lwr)^5)
S7::method(dlinkinv, LowerBoundedLink) <- function(x, eta) exp_floored(eta)
S7::method(d2linkinv, LowerBoundedLink) <- function(x, eta) exp_floored(eta)
S7::method(d3linkinv, LowerBoundedLink) <- function(x, eta) exp_floored(eta)
S7::method(d4linkinv, LowerBoundedLink) <- function(x, eta) exp_floored(eta)
S7::method(d5linkinv, LowerBoundedLink) <- function(x, eta) exp_floored(eta)

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
S7::method(d5linkfun, UpperBoundedLink) <- function(x, theta) -24 / ((x@upr - theta)^5)
S7::method(dlinkinv, UpperBoundedLink) <- function(x, eta) -exp_floored(eta)
S7::method(d2linkinv, UpperBoundedLink) <- function(x, eta) -exp_floored(eta)
S7::method(d3linkinv, UpperBoundedLink) <- function(x, eta) -exp_floored(eta)
S7::method(d4linkinv, UpperBoundedLink) <- function(x, eta) -exp_floored(eta)
S7::method(d5linkinv, UpperBoundedLink) <- function(x, eta) -exp_floored(eta)

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
#' @param lwr Numeric or `NULL`. The lower bound of the interval.
#' @param upr Numeric or `NULL`. The upper bound of the interval.
#'
#' @details
#' The link depends on which endpoints are given:
#'
#' - `lwr` and `upr`: \eqn{\theta} is mapped to
#'   \eqn{p = (\theta - \text{lwr})/(\text{upr} - \text{lwr})} in \eqn{(0, 1)},
#'   and then through the logit;
#' - `lwr` only: \eqn{\eta = \log(\theta - \text{lwr})}, with inverse
#'   \eqn{\theta = \exp(\eta) + \text{lwr}};
#' - `upr` only: \eqn{\eta = \log(\text{upr} - \theta)}, with inverse
#'   \eqn{\theta = \text{upr} - \exp(\eta)};
#' - neither: the result is [identity_link()].
#'
#' @return An S7 object inheriting from [link()], whose methods compute the
#'   link, its inverse and their derivatives to the fifth order. Its class
#'   depends on the endpoints given: [DoublyBoundedLink()],
#'   [LowerBoundedLink()], [UpperBoundedLink()], or [IdentityLink()] when
#'   neither endpoint is supplied.
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
#' @seealso [link()], [logit_link()], [log_link()]
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

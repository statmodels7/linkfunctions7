#' @include generics.R
#' @include link_class.R
NULL

# Numerical fallbacks for the derivative generics.
#
# A link is a formula and its derivatives, so writing all ten methods is a dozen
# lines and the catalogue does it everywhere. But a user experimenting with a
# transformation should not have to differentiate it by hand four times before
# anything works, so the base class supplies every derivative it was not given.
#
# The one rule these obey is the rule stated in check_link(): never differentiate
# numerically more than once in a row. Each fallback finds the highest order the
# link supplies analytically and applies a single stencil of the remaining order
# to *that*, rather than composing four first derivatives. A link that stops at
# the second order therefore gets a third and fourth of ordinary quality, and
# only a link supplying nothing but linkfun() falls back on a fourth-order
# stencil applied to the function itself -- which is honest arithmetic, but
# arithmetic whose accuracy is reported rather than assumed.


#' Highest Analytically Implemented Derivative Order
#'
#' @description
#' The largest \eqn{k \le 4} for which the link's own class registers a method
#' for the order-\eqn{k} derivative generic; \code{0} when it registers none, so
#' that only \code{\link{linkfun}} or \code{\link{linkinv}} is available.
#'
#' @details
#' Detection uses the documented S7 property that a method records the class it
#' was registered on in its \code{signature} attribute: a method inherited from
#' the base \code{\link{link}} class is a fallback, anything else is the link's
#' own. Comparing method objects with \code{identical()} does not work for this,
#' because S7 wraps them.
#'
#' The recorded class is compared by name and package rather than with
#' \code{identical()}, since \code{identical()} on S7 class objects tests
#' object identity and returns \code{FALSE} for a class re-created from the
#' same definition, as happens when the package's code is re-evaluated under
#' coverage instrumentation; identity is kept only as a fast path.
#'
#' The search stops at the first missing order, so the answer always means
#' that every order up to it is analytic.
#'
#' @param x An object of class \code{link}.
#' @param inverse Logical; \code{TRUE} to ask about the inverse-link generics.
#'
#' @return An integer between 0 and 4.
#'
#' @seealso \code{\link{link_fallback_orders}}, which reports this to the user.
#' @keywords internal
analytic_order <- function(x, inverse = FALSE) {
  gens <- if (inverse) {
    list(dlinkinv, d2linkinv, d3linkinv, d4linkinv)
  } else {
    list(dlinkfun, d2linkfun, d3linkfun, d4linkfun)
  }
  cls <- S7::S7_class(x)
  n <- 0L
  for (k in seq_along(gens)) {
    m <- tryCatch(S7::method(gens[[k]], cls), error = function(e) NULL)
    if (is.null(m)) break
    reg <- tryCatch(attr(m, "signature")[[1]], error = function(e) NULL)
    if (is.null(reg) || is_base_link_class(reg)) break
    n <- k
  }
  n
}


#' Is a Class the Base Link Class
#'
#' @description
#' Answers whether an S7 class is this package's own \code{\link{link}} class,
#' which is how a method inherited from the base class is told from one a link
#' registered for itself.
#'
#' @details
#' Object identity is tried first, since it is the usual case and costs
#' nothing, and the class name and package are compared after. The second
#' comparison is what makes the answer reliable: \code{identical()} on an S7
#' class is object identity, so it is false for a class re-created from the
#' same definition, which is what happens whenever a package's code is
#' evaluated rather than loaded. A base fallback mistaken for an analytic
#' method makes every fallback differentiate the order below it, which is the
#' nested differencing the design exists to forbid.
#'
#' @param cls An S7 class.
#'
#' @return \code{TRUE} or \code{FALSE}.
#'
#' @seealso \code{\link{link_fallback_orders}}
#'
#' @keywords internal
is_base_link_class <- function(cls) {
  if (identical(cls, link)) return(TRUE)
  nm <- attr(cls, "name")
  pk <- attr(cls, "package")
  identical(nm, attr(link, "name")) && identical(pk, attr(link, "package"))
}


#' A Finite-Difference Step for a Given Order
#'
#' @description
#' The step for a central stencil of order \code{order}, scaled by the magnitude
#' of the evaluation point and, when bounds are supplied, shrunk so that the
#' whole stencil stays strictly inside them.
#'
#' @details
#' The unclamped step is \eqn{\varepsilon^{1/(k+2)}\max(1, \lvert x\rvert)},
#' which balances truncation against rounding for a central difference of order
#' \eqn{k}: higher orders divide by a higher power of \eqn{h}, so they need a
#' larger one.
#'
#' The clamp matters because a link's domain is open. Differentiating the log
#' link at \eqn{\theta = 10^{-8}} with a step chosen from the magnitude alone
#' would evaluate \eqn{\log} at a negative number; the stencil for order 3 and 4
#' reaches \eqn{2h}, so it is that reach, not \eqn{h}, that must fit.
#'
#' @param x A numeric vector of evaluation points.
#' @param order The derivative order, 1 to 4.
#' @param bounds An optional length-2 numeric vector, the open interval
#'   \code{x} must stay inside.
#'
#' @return A numeric vector of steps, the same length as \code{x}.
#'
#' @keywords internal
fd_step <- function(x, order, bounds = NULL) {
  # The step rule moved to numericals7 with the rest of the stencil library;
  # this wrapper keeps the policy name the fallbacks speak through.
  numericals7::fd_step(x, order, bounds = bounds)
}


#' One Central Stencil, Never Nested
#'
#' @description
#' The order-\code{order} central finite difference of \code{f} at \code{x},
#' applied in a single step rather than by composing lower-order differences.
#'
#' @details
#' The four stencils are
#' \deqn{f' \approx \frac{f(x+h) - f(x-h)}{2h}, \qquad
#'       f'' \approx \frac{f(x+h) - 2f(x) + f(x-h)}{h^{2}},}
#' \deqn{f''' \approx \frac{f(x+2h) - 2f(x+h) + 2f(x-h) - f(x-2h)}{2h^{3}}, \qquad
#'       f'''' \approx \frac{f(x+2h) - 4f(x+h) + 6f(x) - 4f(x-h) + f(x-2h)}{h^{4}}.}
#'
#' Applying one stencil of order \eqn{k} is not the same as applying \eqn{k}
#' stencils of order one, and the difference is the whole reason this function
#' exists: each numerical differentiation multiplies the error of the one before
#' it, so a fourth derivative reached by four nested first differences is noise.
#' The identity link illustrates the failure directly -- its third
#' derivative is exactly zero, and nested differentiation returns a number of
#' order one.
#'
#' @param f A vectorised function of one numeric argument.
#' @param x A numeric vector of evaluation points.
#' @param order The derivative order, 1 to 4.
#' @param h A numeric vector of steps, from \code{\link{fd_step}}.
#'
#' @return A numeric vector of the same length as \code{x}.
#'
#' @keywords internal
stencil_deriv <- function(f, x, order, h) {
  # One stencil of the requested order from numericals7's shared weights;
  # the four formulas the documentation displays are exactly what the
  # Vandermonde construction produces at the default accuracy.
  numericals7::fd_derivative(f, x, order, h = h)
}


#' The Range a Stencil May Evaluate the Inverse Link On
#'
#' @description
#' The image of the link's parameter bounds under \code{\link{linkfun}}, used
#' to keep a finite-difference grid inside the set the inverse link is defined
#' on.
#'
#' @details
#' A link need not map onto the whole real line: the square root reaches only
#' the positive half, and a stencil straying outside returns \code{NaN}, which
#' would make a numerical derivative missing rather than inaccurate. The
#' bounds are returned sorted, since a decreasing link reverses them, and are
#' infinite in the directions where they cannot be established.
#'
#' @param x A \code{\link{link}} object.
#'
#' @return A numeric vector of length two.
#'
#' @seealso \code{\link{link_bounds_clamp}}
#'
#' @keywords internal
eta_bounds <- function(x) {
  b <- tryCatch(sort(linkfun(x, x@link_bounds)), error = function(e) NULL)
  if (is.null(b) || length(b) != 2L || anyNA(b)) c(-Inf, Inf) else b
}


#' The Body Shared by Every Numerical Fallback
#'
#' @description
#' Computes the order-\code{order} derivative of a link, in either direction, by
#' differentiating once the highest order the link supplies analytically.
#'
#' @details
#' The whole design is in the two lines that pick \code{m} and \code{gap}: never
#' differentiate numerically more than the number of orders actually missing.
#' A link analytic to the second order asks for a first difference to reach the
#' third, not three; a link supplying nothing but \code{\link{linkfun}} is the
#' only case in which a fourth-order stencil is applied to the function itself.
#'
#' Note the recursion is only apparent. The base function is fetched through
#' \code{\link{linkderiv}}, which dispatches to the link's own method for an
#' order it implements — so the chain always terminates on analytic code, and
#' never on another fallback.
#'
#' @param x An object of class \code{link}.
#' @param v A numeric vector: \eqn{\theta} going forward, \eqn{\eta} coming back.
#' @param order The derivative order wanted, 1 to 4.
#' @param inverse Logical; \code{TRUE} for the inverse-link direction.
#'
#' @return A numeric vector of the same length as \code{v}.
#'
#' @keywords internal
fallback_deriv <- function(x, v, order, inverse) {
  m <- min(analytic_order(x, inverse = inverse), order - 1L)
  gap <- order - m

  base_fun <- if (inverse) {
    if (m == 0L) function(z) linkinv(x, z) else function(z) linkinvderiv(x, z, order = m)
  } else {
    if (m == 0L) function(z) linkfun(x, z) else function(z) linkderiv(x, z, order = m)
  }

  bnds <- if (inverse) eta_bounds(x) else x@link_bounds
  h <- fd_step(v, gap, bnds)
  na_from(stencil_deriv(base_fun, v, gap, h), v)
}

# S7 requires a method's formals to match the generic's exactly, and the two
# directions name their argument differently, so the wrappers are written out
# rather than generated.
S7::method(dlinkfun,  link) <- function(x, theta) fallback_deriv(x, theta, 1L, FALSE)
S7::method(d2linkfun, link) <- function(x, theta) fallback_deriv(x, theta, 2L, FALSE)
S7::method(d3linkfun, link) <- function(x, theta) fallback_deriv(x, theta, 3L, FALSE)
S7::method(d4linkfun, link) <- function(x, theta) fallback_deriv(x, theta, 4L, FALSE)

S7::method(dlinkinv,  link) <- function(x, eta) fallback_deriv(x, eta, 1L, TRUE)
S7::method(d2linkinv, link) <- function(x, eta) fallback_deriv(x, eta, 2L, TRUE)
S7::method(d3linkinv, link) <- function(x, eta) fallback_deriv(x, eta, 3L, TRUE)
S7::method(d4linkinv, link) <- function(x, eta) fallback_deriv(x, eta, 4L, TRUE)


#' Which Derivative Orders a Link Computes Exactly
#'
#' @description
#' Reports, for each direction, how many derivative orders the link implements
#' analytically and which are therefore obtained by finite differences.
#'
#' @details
#' Every link can answer every derivative generic, because the base class
#' supplies numerical fallbacks for the orders a link does not implement. That
#' convenience requires a way of asking which is which — a fallback is
#' correct but not exact, and it is the reason \code{\link{check_link}} reports
#' such orders separately rather than passing them.
#'
#' @param x An object of class \code{link}.
#'
#' @return A list with \code{forward} and \code{inverse}, each an integer: the
#'   number of leading orders implemented analytically, from 0 to 4.
#'
#' @examples
#' # everything the package ships is exact to fourth order
#' link_fallback_orders(logit_link())
#'
#' @seealso \code{\link{check_link}}
#' @export
link_fallback_orders <- function(x) {
  list(forward = analytic_order(x, inverse = FALSE),
       inverse = analytic_order(x, inverse = TRUE))
}

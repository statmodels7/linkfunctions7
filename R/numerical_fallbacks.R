#' @include generics.R
#' @include link_class.R
NULL

# Numerical fallbacks for the derivative generics.
#
# A link is a formula and its derivatives, so writing all ten methods is a dozen
# lines and the catalog does it everywhere. But a user experimenting with a
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
#' The largest \eqn{k \le 5} for which the link's own class registers a method
#' for the order-\eqn{k} derivative generic; `0` when it registers none, so
#' that only [linkfun()] or [linkinv()] is available.
#'
#' @details
#' Detection uses the documented S7 property that a method records, in its
#' `signature` attribute, the class on which it was registered: a method inherited from
#' the base [link()] class is a fallback, anything else is the link's
#' own. Comparing method objects with `identical()` does not work for this,
#' because S7 wraps them.
#'
#' The recorded class is compared by name and package. `identical()` on S7
#' class objects tests object identity and returns `FALSE` for a class
#' re-created from the same definition, as happens when the package's code is
#' re-evaluated under coverage instrumentation, so identity serves only as a
#' fast path.
#'
#' The search stops at the first missing order, so the answer always means
#' that every order up to it is analytic.
#'
#' @param x An object of class `link`.
#' @param inverse Logical; `TRUE` for the inverse-link generics.
#'
#' @return An integer between 0 and 5.
#'
#' @seealso [link_fallback_orders()], which reports this to the user.
#' @keywords internal
analytic_order <- function(x, inverse = FALSE) {
  gens <- if (inverse) {
    list(dlinkinv, d2linkinv, d3linkinv, d4linkinv, d5linkinv)
  } else {
    list(dlinkfun, d2linkfun, d3linkfun, d4linkfun, d5linkfun)
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
#' Returns whether an S7 class is this package's own [link()] class. It is
#' used to tell a method inherited from the base class from one that a link
#' registered for itself.
#'
#' @details
#' Object identity is tried first, since it is the usual case and costs
#' nothing, and the class name and package are compared after. The second
#' comparison is the one that makes the answer reliable. `identical()` on an
#' S7 class is object identity, so it is false for a class re-created from the
#' same definition, which happens whenever a package's code is evaluated
#' instead of loaded, as it is under coverage instrumentation. If a base
#' fallback were mistaken for an analytic method, every fallback would
#' differentiate the order below it, which is the nested differencing that the
#' fallbacks are written to avoid.
#'
#' @param cls An S7 class.
#'
#' @return `TRUE` or `FALSE`.
#'
#' @seealso [link_fallback_orders()]
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
#' The step for a central stencil of order `order`, scaled by the magnitude
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
#' @param order The derivative order, 1 to 5.
#' @param bounds An optional length-2 numeric vector, the open interval
#'   inside which `x` must stay.
#'
#' @return A numeric vector of steps, the same length as `x`.
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
#' The order-`order` central finite difference of `f` at `x`,
#' applied in a single step rather than by composing lower-order differences.
#'
#' @details
#' The stencils of orders one to four are
#' \deqn{f' \approx \frac{f(x+h) - f(x-h)}{2h}, \qquad
#'       f'' \approx \frac{f(x+h) - 2f(x) + f(x-h)}{h^{2}},}
#' \deqn{f''' \approx \frac{f(x+2h) - 2f(x+h) + 2f(x-h) - f(x-2h)}{2h^{3}}, \qquad
#'       f'''' \approx \frac{f(x+2h) - 4f(x+h) + 6f(x) - 4f(x-h) + f(x-2h)}{h^{4}}.}
#'
#' Applying one stencil of order \eqn{k} differs from applying \eqn{k}
#' stencils of order one, and this function exists because of the difference:
#' each numerical differentiation amplifies the error of the one before it, so
#' a fourth derivative computed by four nested first differences is dominated
#' by rounding error. For the identity link, whose third derivative is exactly
#' zero, nested differentiation returns non-zero values.
#'
#' @param f A vectorized function of one numeric argument.
#' @param x A numeric vector of evaluation points.
#' @param order The derivative order, 1 to 5.
#' @param h A numeric vector of steps, from [fd_step()].
#'
#' @return A numeric vector of the same length as `x`.
#'
#' @keywords internal
stencil_deriv <- function(f, x, order, h) {
  # One stencil of the requested order from numericals7's shared weights;
  # the four formulas the documentation displays are exactly what the
  # Vandermonde construction produces at the default accuracy.
  numericals7::fd_derivative(f, x, order, h = h)
}


#' The Range of Predictors a Link Admits
#'
#' @description
#' The image of the link's parameter bounds under [linkfun()], which is the
#' set of predictors on which the inverse link is defined. Where a link maps
#' onto the whole real line the result is `c(-Inf, Inf)`; otherwise its finite
#' end or ends bound the values that a caller may pass to [linkinv()].
#'
#' @details
#' A link need not map onto the whole real line: the square root reaches only
#' the positive half, and so do [inverse_link()], [inverse_sq_link()] and
#' [power_link()] at any non-zero exponent. The bounds are returned sorted, since
#' a decreasing link reverses them, and are infinite in the directions where
#' they cannot be established.
#'
#' # Use outside the package
#'
#' Inside the package the function keeps a finite-difference grid within the
#' set on which the inverse link is defined, since a stencil point outside it
#' returns `NaN` and makes a numerical derivative missing. Outside the package
#' it serves a related purpose: [link_bounds()][link] gives the set that a link
#' maps **onto**, and a consumer that carries an unconstrained vector needs the
#' set that it maps **from**.
#'
#' A family that reads a free vector in \eqn{\mathbb{R}^d} and applies an
#' inverse link to each coordinate needs the map to be defined and injective
#' on the whole real line, which fails for a link with finite eta bounds:
#' [linkinv()] of a square root link is even, so `-2` and `2` both give 4 and
#' the round trip returns the absolute value. Such a family can reject the
#' link at construction by testing `all(is.infinite(eta_bounds(link)))`, and
#' its error message can then name the link.
#'
#' @param x A [link()] object.
#'
#' @return A numeric vector of length two, sorted, with `-Inf` or `Inf` in
#'   either position where that end is unbounded.
#'
#' @seealso [link_bounds_clamp()], which keeps [linkinv()]'s result strictly
#'   inside the parameter bounds, and [linkfun()], the map whose image this is.
#'
#' @examples
#' # A link from the whole real line onto the positive half of the theta axis.
#' eta_bounds(log_link())
#'
#' # The square root reaches only the positive half of the eta axis, so its
#' # inverse is even there and the round trip returns the absolute value.
#' eta_bounds(sqrt_link())
#' linkfun(sqrt_link(), linkinv(sqrt_link(), -2))
#'
#' # Which is the question a consumer carrying an unconstrained vector asks.
#' links <- list(log_link(), softplus_link(), logit_link(), sqrt_link(),
#'               inverse_link(), inverse_sq_link(), power_link(0.5))
#' vapply(links, function(l) all(is.infinite(eta_bounds(l))), logical(1))
#'
#' @export
eta_bounds <- function(x) {
  b <- tryCatch(sort(linkfun(x, x@link_bounds)), error = function(e) NULL)
  if (is.null(b) || length(b) != 2L || anyNA(b)) c(-Inf, Inf) else b
}


#' The Body Shared by Every Numerical Fallback
#'
#' @description
#' Computes the order-`order` derivative of a link, in either direction, by
#' one numerical differentiation of the highest order that the link supplies
#' analytically.
#'
#' @details
#' The function takes the highest analytic order \eqn{m} below the requested
#' one and applies a single central stencil of order `order - m` to it, never
#' a chain of lower-order stencils. For a link analytic to the second order,
#' the third derivative is one first difference of the second; only for a link
#' that supplies nothing but [linkfun()] is a stencil of the full order
#' applied to the function itself.
#'
#' The recursion is only apparent. The base function is obtained through
#' [linkderiv()] or [linkinvderiv()], which dispatch to the link's own method
#' for an order that
#' the link implements, so the chain always ends on analytic code and never on
#' another fallback.
#'
#' @param x An object of class `link`.
#' @param v A numeric vector: \eqn{\theta} going forward, \eqn{\eta} coming back.
#' @param order The derivative order wanted, 1 to 5.
#' @param inverse Logical; `TRUE` for the inverse-link direction.
#'
#' @return A numeric vector of the same length as `v`.
#'
#' @section Methods:
#' The ten registrations on the base class [link()] have this function as
#' their whole body: `dlinkfun()` through `d5linkfun()` going out and
#' `dlinkinv()` through `d5linkinv()` coming back, each passing its order and
#' its direction. A link inherits them for the orders it does not implement
#' itself, so a link defined with nothing but [linkfun()] and [linkinv()]
#' still has a method for every derivative generic. S7 requires the formals of
#' a method to match those of the generic, and the two directions name their
#' argument differently, so the ten wrappers are written out instead of
#' generated.
#'
#' @aliases dlinkfun.link
#' @aliases d2linkfun.link
#' @aliases d3linkfun.link
#' @aliases d4linkfun.link
#' @aliases d5linkfun.link
#' @aliases dlinkinv.link
#' @aliases d2linkinv.link
#' @aliases d3linkinv.link
#' @aliases d4linkinv.link
#' @aliases d5linkinv.link
#'
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
S7::method(d5linkfun, link) <- function(x, theta) fallback_deriv(x, theta, 5L, FALSE)

S7::method(dlinkinv,  link) <- function(x, eta) fallback_deriv(x, eta, 1L, TRUE)
S7::method(d2linkinv, link) <- function(x, eta) fallback_deriv(x, eta, 2L, TRUE)
S7::method(d3linkinv, link) <- function(x, eta) fallback_deriv(x, eta, 3L, TRUE)
S7::method(d4linkinv, link) <- function(x, eta) fallback_deriv(x, eta, 4L, TRUE)
S7::method(d5linkinv, link) <- function(x, eta) fallback_deriv(x, eta, 5L, TRUE)


#' Which Derivative Orders a Link Computes Exactly
#'
#' @description
#' Reports, for each direction, how many derivative orders the link implements
#' analytically and which are therefore obtained by finite differences.
#'
#' @details
#' Every link has a method for every derivative generic, because the base
#' class supplies numerical fallbacks for the orders that a link does not
#' implement. This function reports which orders are analytic. A fallback is
#' accurate only to the precision of a finite difference, so [check_link()]
#' reports such orders separately instead of counting them as passed.
#'
#' @param x An object of class `link`.
#'
#' @return A list with `forward` and `inverse`, each an integer: the
#'   number of leading orders implemented analytically, from 0 to 5.
#'
#' @examples
#' # everything the package ships is exact to fifth order
#' link_fallback_orders(logit_link())
#'
#' @seealso [check_link()]
#' @export
link_fallback_orders <- function(x) {
  list(forward = analytic_order(x, inverse = FALSE),
       inverse = analytic_order(x, inverse = TRUE))
}

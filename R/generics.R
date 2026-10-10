#' @title Evaluate Forward Link Function
#' @include link_class.R
#'
#' @description
#' Carries a parameter from its own domain onto the real line,
#' \eqn{\eta = g(\theta)}.
#'
#' @details
#' A link is a strictly monotone differentiable bijection
#' \eqn{g : \Theta \to \mathbb{R}} from the parameter's open domain
#' `x@link_bounds` onto the whole line, so that an unconstrained
#' optimizer may work in \eqn{\eta} while \eqn{\theta = g^{-1}(\eta)} stays
#' admissible at every point. [linkinv()] evaluates that inverse
#' and is the only other method that a link must supply; the derivative
#' generics have numerical fallbacks derived from the pair.
#'
#' @param x An object of class `link`.
#' @param theta A numeric vector of parameters.
#' @return A numeric vector of the linear predictor.
#' @examples
#' linkfun(logit_link(), c(0.25, 0.5, 0.75))
#' linkfun(log_link(), c(1, exp(1)))
#' @seealso [linkinv()], [linkderiv()], [linkinvderiv()], [check_link()], [link()]
#' @export
linkfun <- S7::new_generic("linkfun", "x", fun = function(x, theta) S7::S7_dispatch())

#' @title Evaluate Inverse Link Function
#'
#' @description
#' Maps a linear predictor back onto the parameter's domain.
#'
#' @details
#' A link is a bijection onto an open interval, so that a value returned by
#' the inverse link can be passed back to [linkfun()], or to a density that
#' validates its parameters against open intervals. In double precision the
#' map is not quite onto: `plogis` is exactly 1 above about \eqn{\eta = 37},
#' `lwr + exp(eta)` rounds to `lwr` once the exponential falls below half an
#' ulp of a non-zero `lwr`, both overflow to infinity for large arguments, and
#' the square-root, power, inverse and inverse-square links reach a bound or
#' overflow at \eqn{\eta = 0}.
#'
#' The generic therefore clamps the result with [link_bounds_clamp()]: a value
#' exactly equal to a finite bound is moved one or two units in the last place
#' inside it (to the smallest positive normal double when the bound is zero),
#' and an infinite value becomes the largest finite double of that sign. A
#' value outside the bounds by more than rounding, which comes from a linear
#' predictor outside [eta_bounds()], is returned unchanged, and so are `NA`
#' and `NaN`. The clamp costs one comparison per bound.
#'
#' Because the clamp is applied in the generic, every link inherits it,
#' including a user-defined one, and a method can be written as the plain
#' formula without a guard of its own.
#'
#' @param x An object of class `link`.
#' @param eta A numeric vector of linear predictors.
#'
#' @return A numeric vector of the same length as `eta`. A value that lands on
#'   a bound or overflows is moved strictly inside `x@link_bounds`; a value
#'   outside the bounds by more than rounding, `NA` and `NaN` are returned
#'   unchanged.
#'
#' @examples
#' linkinv(logit_link(), c(-1, 0, 1))
#' linkinv(log_link(), c(0, 1))
#'
#' # far out, where the arithmetic saturates: strictly inside (0, 1) rather
#' # than exactly 1, so the round trip still returns a number
#' linkinv(logit_link(), 40) < 1
#' linkfun(logit_link(), linkinv(logit_link(), 40))
#'
#' @seealso [linkfun()], [link_bounds_clamp()]
#' @export
linkinv <- S7::new_generic("linkinv", "x", fun = function(x, eta) {
  # S7_dispatch() has to be assigned on its own line rather than nested in the
  # call: evaluated as an argument it re-enters the promise for `eta` and fails
  # with "recursive default argument reference". The same shape as
  # distributions7's link-scale interception, and for the same reason.
  theta <- S7::S7_dispatch()
  link_bounds_clamp(theta, x@link_bounds)
})

#' @title 1st Derivative of a Link Function
#' @description
#' The first derivative of the link \eqn{g(\theta)} with respect to the
#' parameter, on the parameter scale. A delta-method standard error is carried
#' in this direction, from a variance on \eqn{\theta} to one on \eqn{\eta}.
#' @details
#' Every link has a method for this generic. If the class of a link registers
#' no method for it, the numerical method of the base class is used: it
#' applies a single central stencil to the highest order that the link
#' supplies analytically, never a chain of lower-order differences.
#' [link_fallback_orders()] reports which orders of a given link are exact,
#' and [check_link()] leaves a fallback order unchecked, because comparing it
#' with a finite difference of the order below would repeat the same
#' computation and could not detect an error.
#'
#' In code where speed matters, call this generic directly: [linkderiv()] and
#' [linkinvderiv()] route by order and so dispatch twice, once on themselves
#' and once here, and on a short vector the routed call takes two to three
#' times as long as the direct one.
#' @param x An object of class `link`.
#' @param theta A numeric vector of parameter values, inside
#'   `x@link_bounds`. The domain is not checked: outside it the formula is
#'   evaluated as written, so the result may be `NaN`, `NA` or an ordinary
#'   number, and no error is signaled.
#' @return A numeric vector of the same length as `theta`, missing wherever
#'   `theta` is.
#' @seealso [linkderiv()], which routes to this generic by order, and
#'   [dlinkinv()] for the same order in the other direction.
#' @examples
#' # The log link's forward derivatives are 1 / t, so at theta = 2:
#' dlinkfun(log_link(), 2) - (1 / 2)
#'
#' # Missingness propagates instead of being filled in.
#' dlinkfun(logit_link(), c(0.5, NA))
#' @export
dlinkfun <- S7::new_generic("dlinkfun", "x", fun = function(x, theta) S7::S7_dispatch())

#' @title 2nd Derivative of a Link Function
#' @description
#' The second derivative of the link \eqn{g(\theta)} with respect to the
#' parameter, on the parameter scale.
#' @details
#' Every link has a method for this generic. If the class of a link registers
#' no method for it, the numerical method of the base class is used: it
#' applies a single central stencil to the highest order that the link
#' supplies analytically, never a chain of lower-order differences.
#' [link_fallback_orders()] reports which orders of a given link are exact,
#' and [check_link()] leaves a fallback order unchecked, because comparing it
#' with a finite difference of the order below would repeat the same
#' computation and could not detect an error.
#'
#' In code where speed matters, call this generic directly: [linkderiv()] and
#' [linkinvderiv()] route by order and so dispatch twice, once on themselves
#' and once here, and on a short vector the routed call takes two to three
#' times as long as the direct one.
#' @param x An object of class `link`.
#' @param theta A numeric vector of parameter values, inside
#'   `x@link_bounds`. The domain is not checked: outside it the formula is
#'   evaluated as written, so the result may be `NaN`, `NA` or an ordinary
#'   number, and no error is signaled.
#' @return A numeric vector of the same length as `theta`, missing wherever
#'   `theta` is.
#' @seealso [linkderiv()], which routes to this generic by order, and
#'   [d2linkinv()] for the same order in the other direction.
#' @examples
#' # The log link's forward derivatives are -1 / t^2, so at theta = 2:
#' d2linkfun(log_link(), 2) - (-1 / 2^2)
#'
#' # Missingness propagates instead of being filled in.
#' d2linkfun(logit_link(), c(0.5, NA))
#' @export
d2linkfun <- S7::new_generic("d2linkfun", "x", fun = function(x, theta) S7::S7_dispatch())

#' @title 3rd Derivative of a Link Function
#' @description
#' The third derivative of the link \eqn{g(\theta)} with respect to the
#' parameter, on the parameter scale.
#' @details
#' Every link has a method for this generic. If the class of a link registers
#' no method for it, the numerical method of the base class is used: it
#' applies a single central stencil to the highest order that the link
#' supplies analytically, never a chain of lower-order differences.
#' [link_fallback_orders()] reports which orders of a given link are exact,
#' and [check_link()] leaves a fallback order unchecked, because comparing it
#' with a finite difference of the order below would repeat the same
#' computation and could not detect an error.
#'
#' In code where speed matters, call this generic directly: [linkderiv()] and
#' [linkinvderiv()] route by order and so dispatch twice, once on themselves
#' and once here, and on a short vector the routed call takes two to three
#' times as long as the direct one.
#' @param x An object of class `link`.
#' @param theta A numeric vector of parameter values, inside
#'   `x@link_bounds`. The domain is not checked: outside it the formula is
#'   evaluated as written, so the result may be `NaN`, `NA` or an ordinary
#'   number, and no error is signaled.
#' @return A numeric vector of the same length as `theta`, missing wherever
#'   `theta` is.
#' @seealso [linkderiv()], which routes to this generic by order, and
#'   [d3linkinv()] for the same order in the other direction.
#' @examples
#' # The log link's forward derivatives are 2 / t^3, so at theta = 2:
#' d3linkfun(log_link(), 2) - (2 / 2^3)
#'
#' # Missingness propagates instead of being filled in.
#' d3linkfun(logit_link(), c(0.5, NA))
#' @export
d3linkfun <- S7::new_generic("d3linkfun", "x", fun = function(x, theta) S7::S7_dispatch())

#' @title 4th Derivative of a Link Function
#' @description
#' The fourth derivative of the link \eqn{g(\theta)} with respect to the
#' parameter, on the parameter scale.
#' @details
#' Every link has a method for this generic. If the class of a link registers
#' no method for it, the numerical method of the base class is used: it
#' applies a single central stencil to the highest order that the link
#' supplies analytically, never a chain of lower-order differences.
#' [link_fallback_orders()] reports which orders of a given link are exact,
#' and [check_link()] leaves a fallback order unchecked, because comparing it
#' with a finite difference of the order below would repeat the same
#' computation and could not detect an error.
#'
#' In code where speed matters, call this generic directly: [linkderiv()] and
#' [linkinvderiv()] route by order and so dispatch twice, once on themselves
#' and once here, and on a short vector the routed call takes two to three
#' times as long as the direct one.
#' @param x An object of class `link`.
#' @param theta A numeric vector of parameter values, inside
#'   `x@link_bounds`. The domain is not checked: outside it the formula is
#'   evaluated as written, so the result may be `NaN`, `NA` or an ordinary
#'   number, and no error is signaled.
#' @return A numeric vector of the same length as `theta`, missing wherever
#'   `theta` is.
#' @seealso [linkderiv()], which routes to this generic by order, and
#'   [d4linkinv()] for the same order in the other direction.
#' @examples
#' # The log link's forward derivatives are -6 / t^4, so at theta = 2:
#' d4linkfun(log_link(), 2) - (-6 / 2^4)
#'
#' # Missingness propagates instead of being filled in.
#' d4linkfun(logit_link(), c(0.5, NA))
#' @export
d4linkfun <- S7::new_generic("d4linkfun", "x", fun = function(x, theta) S7::S7_dispatch())
#' @title 5th Derivative of a Link Function
#' @description
#' The fifth derivative of the link \eqn{g(\theta)} with respect to the
#' parameter, on the parameter scale.
#' @details
#' The fifth order is used by score-driven filters. Each order of
#' differentiation of a filtered predictor through its recursion involves one
#' more order of the link and of the family, because the score that drives the
#' recursion is evaluated at the predictor that the recursion produces.
#'
#' Every link has a method for this generic. If the class of a link registers
#' no method for it, the numerical method of the base class is used: it
#' applies a single central stencil to the highest order that the link
#' supplies analytically, never a chain of lower-order differences.
#' [link_fallback_orders()] reports which orders of a given link are exact,
#' and [check_link()] leaves a fallback order unchecked, because comparing it
#' with a finite difference of the order below would repeat the same
#' computation and could not detect an error.
#'
#' In code where speed matters, call this generic directly: [linkderiv()] and
#' [linkinvderiv()] route by order and so dispatch twice, once on themselves
#' and once here, and on a short vector the routed call takes two to three
#' times as long as the direct one.
#' @param x An object of class `link`.
#' @param theta A numeric vector of parameter values, inside
#'   `x@link_bounds`. The domain is not checked: outside it the formula is
#'   evaluated as written, so the result may be `NaN`, `NA` or an ordinary
#'   number, and no error is signaled.
#' @return A numeric vector of the same length as `theta`, missing wherever
#'   `theta` is.
#' @seealso [linkderiv()], which routes to this generic by order, and
#'   [d5linkinv()] for the same order in the other direction.
#' @examples
#' # The log link's forward derivatives are 24 / t^5, so at theta = 2:
#' d5linkfun(log_link(), 2) - (24 / 2^5)
#'
#' # Missingness propagates instead of being filled in.
#' d5linkfun(logit_link(), c(0.5, NA))
#' @export
d5linkfun <- S7::new_generic("d5linkfun", "x", fun = function(x, theta) S7::S7_dispatch())

#' @title 1st Derivative of an Inverse Link Function
#' @description
#' The first derivative of the inverse link \eqn{g^{-1}(\eta)} with respect to
#' the linear predictor. A modeling routine that works on the unconstrained
#' scale uses it as a chain-rule factor, which carries a derivative of the
#' log-likelihood from \eqn{\theta} onto \eqn{\eta}.
#' @details
#' Every link has a method for this generic. If the class of a link registers
#' no method for it, the numerical method of the base class is used: it
#' applies a single central stencil to the highest order that the link
#' supplies analytically, never a chain of lower-order differences.
#' [link_fallback_orders()] reports which orders of a given link are exact,
#' and [check_link()] leaves a fallback order unchecked, because comparing it
#' with a finite difference of the order below would repeat the same
#' computation and could not detect an error.
#'
#' In code where speed matters, call this generic directly: [linkderiv()] and
#' [linkinvderiv()] route by order and so dispatch twice, once on themselves
#' and once here, and on a short vector the routed call takes two to three
#' times as long as the direct one.
#' @param x An object of class `link`.
#' @param eta A numeric vector of linear predictors, inside the range that
#'   [eta_bounds()] returns for the link. The derivative is computed from
#'   `eta` directly, without the clamp that [linkinv()] applies.
#' @return A numeric vector of the same length as `eta`, missing wherever
#'   `eta` is.
#' @seealso [linkinvderiv()], which routes to this generic by order, and
#'   [dlinkfun()] for the same order in the other direction.
#' @examples
#' # Every derivative of exp is exp, so the log link's inverse gives the
#' # same number at every order.
#' dlinkinv(log_link(), 1) - exp(1)
#'
#' # The logit's inverse derivatives are polynomials in theta. At eta = 0 the
#' # logistic is symmetric about 1/2, so its even-order derivatives vanish
#' # there while the first is the Bernoulli variance, 1/4.
#' dlinkinv(logit_link(), 0)
#' @export
dlinkinv <- S7::new_generic("dlinkinv", "x", fun = function(x, eta) S7::S7_dispatch())

#' @title 2nd Derivative of an Inverse Link Function
#' @description
#' The second derivative of the inverse link \eqn{g^{-1}(\eta)} with respect to
#' the linear predictor. It enters the chain rule (the formula of Faa di
#' Bruno) that carries the higher derivatives of the log-likelihood from
#' \eqn{\theta} onto \eqn{\eta}.
#' @details
#' Every link has a method for this generic. If the class of a link registers
#' no method for it, the numerical method of the base class is used: it
#' applies a single central stencil to the highest order that the link
#' supplies analytically, never a chain of lower-order differences.
#' [link_fallback_orders()] reports which orders of a given link are exact,
#' and [check_link()] leaves a fallback order unchecked, because comparing it
#' with a finite difference of the order below would repeat the same
#' computation and could not detect an error.
#'
#' In code where speed matters, call this generic directly: [linkderiv()] and
#' [linkinvderiv()] route by order and so dispatch twice, once on themselves
#' and once here, and on a short vector the routed call takes two to three
#' times as long as the direct one.
#' @param x An object of class `link`.
#' @param eta A numeric vector of linear predictors, inside the range that
#'   [eta_bounds()] returns for the link. The derivative is computed from
#'   `eta` directly, without the clamp that [linkinv()] applies.
#' @return A numeric vector of the same length as `eta`, missing wherever
#'   `eta` is.
#' @seealso [linkinvderiv()], which routes to this generic by order, and
#'   [d2linkfun()] for the same order in the other direction.
#' @examples
#' # Every derivative of exp is exp, so the log link's inverse gives the
#' # same number at every order.
#' d2linkinv(log_link(), 1) - exp(1)
#'
#' # The logit's inverse derivatives are polynomials in theta. At eta = 0 the
#' # logistic is symmetric about 1/2, so its even-order derivatives vanish
#' # there while the first is the Bernoulli variance, 1/4.
#' d2linkinv(logit_link(), 0)
#' @export
d2linkinv <- S7::new_generic("d2linkinv", "x", fun = function(x, eta) S7::S7_dispatch())

#' @title 3rd Derivative of an Inverse Link Function
#' @description
#' The third derivative of the inverse link \eqn{g^{-1}(\eta)} with respect to
#' the linear predictor. It enters the chain rule (the formula of Faa di
#' Bruno) that carries the higher derivatives of the log-likelihood from
#' \eqn{\theta} onto \eqn{\eta}.
#' @details
#' Every link has a method for this generic. If the class of a link registers
#' no method for it, the numerical method of the base class is used: it
#' applies a single central stencil to the highest order that the link
#' supplies analytically, never a chain of lower-order differences.
#' [link_fallback_orders()] reports which orders of a given link are exact,
#' and [check_link()] leaves a fallback order unchecked, because comparing it
#' with a finite difference of the order below would repeat the same
#' computation and could not detect an error.
#'
#' In code where speed matters, call this generic directly: [linkderiv()] and
#' [linkinvderiv()] route by order and so dispatch twice, once on themselves
#' and once here, and on a short vector the routed call takes two to three
#' times as long as the direct one.
#' @param x An object of class `link`.
#' @param eta A numeric vector of linear predictors, inside the range that
#'   [eta_bounds()] returns for the link. The derivative is computed from
#'   `eta` directly, without the clamp that [linkinv()] applies.
#' @return A numeric vector of the same length as `eta`, missing wherever
#'   `eta` is.
#' @seealso [linkinvderiv()], which routes to this generic by order, and
#'   [d3linkfun()] for the same order in the other direction.
#' @examples
#' # Every derivative of exp is exp, so the log link's inverse gives the
#' # same number at every order.
#' d3linkinv(log_link(), 1) - exp(1)
#'
#' # The logit's inverse derivatives are polynomials in theta. At eta = 0 the
#' # logistic is symmetric about 1/2, so its even-order derivatives vanish
#' # there while the first is the Bernoulli variance, 1/4.
#' d3linkinv(logit_link(), 0)
#' @export
d3linkinv <- S7::new_generic("d3linkinv", "x", fun = function(x, eta) S7::S7_dispatch())

#' @title 4th Derivative of an Inverse Link Function
#' @description
#' The fourth derivative of the inverse link \eqn{g^{-1}(\eta)} with respect to
#' the linear predictor. It enters the chain rule (the formula of Faa di
#' Bruno) that carries the higher derivatives of the log-likelihood from
#' \eqn{\theta} onto \eqn{\eta}.
#' @details
#' Every link has a method for this generic. If the class of a link registers
#' no method for it, the numerical method of the base class is used: it
#' applies a single central stencil to the highest order that the link
#' supplies analytically, never a chain of lower-order differences.
#' [link_fallback_orders()] reports which orders of a given link are exact,
#' and [check_link()] leaves a fallback order unchecked, because comparing it
#' with a finite difference of the order below would repeat the same
#' computation and could not detect an error.
#'
#' In code where speed matters, call this generic directly: [linkderiv()] and
#' [linkinvderiv()] route by order and so dispatch twice, once on themselves
#' and once here, and on a short vector the routed call takes two to three
#' times as long as the direct one.
#' @param x An object of class `link`.
#' @param eta A numeric vector of linear predictors, inside the range that
#'   [eta_bounds()] returns for the link. The derivative is computed from
#'   `eta` directly, without the clamp that [linkinv()] applies.
#' @return A numeric vector of the same length as `eta`, missing wherever
#'   `eta` is.
#' @seealso [linkinvderiv()], which routes to this generic by order, and
#'   [d4linkfun()] for the same order in the other direction.
#' @examples
#' # Every derivative of exp is exp, so the log link's inverse gives the
#' # same number at every order.
#' d4linkinv(log_link(), 1) - exp(1)
#'
#' # The logit's inverse derivatives are polynomials in theta. At eta = 0 the
#' # logistic is symmetric about 1/2, so its even-order derivatives vanish
#' # there while the first is the Bernoulli variance, 1/4.
#' d4linkinv(logit_link(), 0)
#' @export
d4linkinv <- S7::new_generic("d4linkinv", "x", fun = function(x, eta) S7::S7_dispatch())
#' @title 5th Derivative of an Inverse Link Function
#' @description
#' The fifth derivative of the inverse link \eqn{g^{-1}(\eta)} with respect
#' to the linear predictor.
#' @details
#' This order is reached by the outer curvature of a score-driven filter. The
#' chain rule that carries a log-likelihood derivative of order \eqn{k} from
#' \eqn{\theta} onto \eqn{\eta} combines the derivatives of the family with
#' the partial Bell polynomials in \eqn{h', \ldots, h^{(k)}}, so an order-5
#' quantity on the unconstrained scale requires the fifth derivative of the
#' inverse link and no higher one.
#'
#' Every link has a method for this generic. If the class of a link registers
#' no method for it, the numerical method of the base class is used: it
#' applies a single central stencil to the highest order that the link
#' supplies analytically, never a chain of lower-order differences.
#' [link_fallback_orders()] reports which orders of a given link are exact,
#' and [check_link()] leaves a fallback order unchecked, because comparing it
#' with a finite difference of the order below would repeat the same
#' computation and could not detect an error.
#'
#' In code where speed matters, call this generic directly: [linkderiv()] and
#' [linkinvderiv()] route by order and so dispatch twice, once on themselves
#' and once here, and on a short vector the routed call takes two to three
#' times as long as the direct one.
#' @param x An object of class `link`.
#' @param eta A numeric vector of linear predictors, inside the range that
#'   [eta_bounds()] returns for the link. The derivative is computed from
#'   `eta` directly, without the clamp that [linkinv()] applies.
#' @return A numeric vector of the same length as `eta`, missing wherever
#'   `eta` is.
#' @seealso [linkinvderiv()], which routes to this generic by order, and
#'   [d5linkfun()] for the same order in the other direction.
#' @examples
#' # Every derivative of exp is exp, so the log link's inverse gives the
#' # same number at every order.
#' d5linkinv(log_link(), 1) - exp(1)
#'
#' # At eta = 0 the logistic is symmetric about 1/2, so its even-order
#' # derivatives vanish there while the odd ones do not.
#' d5linkinv(logit_link(), 0)
#' @export
d5linkinv <- S7::new_generic("d5linkinv", "x", fun = function(x, eta) S7::S7_dispatch())

#' @title Evaluate Derivative of Link Function by Order
#' @description
#' A convenience router over [linkfun()] and the five
#' `d*linkfun` generics.
#' @param x An object of class `link`.
#' @param theta A numeric vector.
#' @param order An integer specifying the derivative order (0 to 5). Order 0 is
#'   the link function itself.
#' @return A numeric vector of the same length as `theta`: the requested
#'   derivative of \eqn{g} evaluated there.
#' @details
#' This dispatches twice, once on itself and once on the order-specific generic.
#' Where that matters, call [dlinkfun()] and its siblings directly.
#' @examples
#' lk <- logit_link()
#' linkderiv(lk, 0.5, order = 0)   # the link itself
#' linkderiv(lk, 0.5, order = 1)
#'
#' # every order at once
#' vapply(0:4, function(k) linkderiv(lk, 0.25, order = k), numeric(1))
#' @seealso [linkfun()], [linkinv()], [linkinvderiv()], [check_link()], [link()]
#' @export
linkderiv <- S7::new_generic("linkderiv", "x", fun = function(x, theta, order = 1) S7::S7_dispatch())

#' @title Evaluate Derivative of Inverse Link Function by Order
#' @description
#' A convenience router over [linkinv()] and the five
#' `d*linkinv` generics.
#' @param x An object of class `link`.
#' @param eta A numeric vector.
#' @param order An integer specifying the derivative order (0 to 5). Order 0 is
#'   the inverse link itself.
#' @return A numeric vector of the same length as `eta`: the requested
#'   derivative of \eqn{g^{-1}} evaluated there.
#' @details
#' This dispatches twice, once on itself and once on the order-specific generic.
#' Where that matters, call [dlinkinv()] and its siblings directly.
#' @examples
#' lk <- logit_link()
#' linkinvderiv(lk, 0, order = 0)  # the inverse link itself
#' linkinvderiv(lk, 0, order = 1)
#'
#' vapply(0:4, function(k) linkinvderiv(lk, 0.5, order = k), numeric(1))
#' @seealso [linkfun()], [linkinv()], [linkderiv()], [check_link()], [link()]
#' @export
linkinvderiv <- S7::new_generic("linkinvderiv", "x", fun = function(x, eta, order = 1) S7::S7_dispatch())

#' @title Validate and Check a Link Object
#' @param x An object of class `link`.
#' @param tolerance Numeric tolerance for floating-point comparisons.
#' @param ... Additional arguments passed to methods.
#' @return Invisibly, a named list of check results; see
#'   [check_link.link()] for its shape. Called mainly for the summary
#'   printed to the console.
#' @examples
#' check_link(sqrt_link())
#' @seealso [linkfun()], [linkinv()], [linkderiv()], [linkinvderiv()], [link()]
#' @export
check_link <- S7::new_generic("check_link", "x", fun = function(x, tolerance = 1e-5, ...) S7::S7_dispatch())

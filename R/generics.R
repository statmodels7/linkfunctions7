#' @title Evaluate Forward Link Function
#' @include link_class.R
#' @param x An object of class \code{link}.
#' @param theta A numeric vector of parameters.
#' @return A numeric vector of the linear predictor.
#' @examples
#' linkfun(logit_link(), c(0.25, 0.5, 0.75))
#' linkfun(log_link(), c(1, exp(1)))
#' @export
linkfun <- S7::new_generic("linkfun", "x", fun = function(x, theta) S7::S7_dispatch())

#' @title Evaluate Inverse Link Function
#'
#' @description
#' Maps a linear predictor back onto the parameter's domain.
#'
#' @details
#' \strong{The result is always strictly inside \code{link_bounds}}, and that is
#' enforced here rather than left to each method. A link is a bijection onto an
#' \emph{open} interval, which is exactly what makes it useful — a value it
#' returns can be handed to \code{\link{linkfun}} and come back, or to a density
#' that validates its parameters against open intervals. In double precision the
#' bijection is not quite onto: \code{plogis} is exactly 1 above about
#' \eqn{\eta = 37}, \code{lwr + exp(eta)} rounds to \code{lwr} once the
#' exponential falls below half an ulp of a non-zero \code{lwr}, and both
#' overflow to infinity far out. Nine of the links shipped here reach a bound
#' somewhere in \eqn{\lvert \eta \rvert \le 800}.
#'
#' So the generic clamps: a result at or beyond a finite bound becomes the
#' nearest double strictly inside it, and a non-finite result becomes the
#' largest finite double of that sign. Both are derived rather than chosen —
#' they are the extremes of what "strictly inside, and a usable number" permits,
#' and no tolerance is invented. The clamp costs one comparison per bound and
#' fires only in the tails.
#'
#' Doing it in the generic body means every link inherits it, including one you
#' write yourself, and that a method can be written as the plain mathematical
#' formula without a guard of its own.
#'
#' @param x An object of class \code{link}.
#' @param eta A numeric vector of linear predictors.
#'
#' @return A numeric vector, strictly inside \code{x@link_bounds}.
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
#' @seealso \code{\link{linkfun}}, \code{\link{link_bounds_clamp}}
#' @export
linkinv <- S7::new_generic("linkinv", "x", fun = function(x, eta) {
  # S7_dispatch() has to be assigned on its own line rather than nested in the
  # call: evaluated as an argument it re-enters the promise for `eta` and fails
  # with "recursive default argument reference". The same shape as
  # distributions7's link-scale interception, and for the same reason.
  theta <- S7::S7_dispatch()
  link_bounds_clamp(theta, x@link_bounds)
})

#' @title 1st Derivative of Link Function
#' @description
#' The first derivative of \eqn{g(\theta)} with respect to \eqn{\theta}.
#' @param x An object of class \code{link}.
#' @param theta A numeric vector.
#' @return A numeric vector of the same length as \code{theta}, missing wherever
#'   \code{theta} is.
#' @seealso \code{\link{linkderiv}}, which routes to this generic by order.
#' @examples
#' dlinkfun(logit_link(), 0.5)
#' dlinkfun(log_link(), c(1, 2))
#' @export
dlinkfun <- S7::new_generic("dlinkfun", "x", fun = function(x, theta) S7::S7_dispatch())

#' @title 2nd Derivative of Link Function
#' @description
#' The second derivative of \eqn{g(\theta)} with respect to \eqn{\theta}.
#' @param x An object of class \code{link}.
#' @param theta A numeric vector.
#' @return A numeric vector of the same length as \code{theta}, missing wherever
#'   \code{theta} is.
#' @seealso \code{\link{linkderiv}}, which routes to this generic by order.
#' @examples
#' d2linkfun(logit_link(), 0.5)
#' d2linkfun(log_link(), c(1, 2))
#' @export
d2linkfun <- S7::new_generic("d2linkfun", "x", fun = function(x, theta) S7::S7_dispatch())

#' @title 3rd Derivative of Link Function
#' @description
#' The third derivative of \eqn{g(\theta)} with respect to \eqn{\theta}.
#' @param x An object of class \code{link}.
#' @param theta A numeric vector.
#' @return A numeric vector of the same length as \code{theta}, missing wherever
#'   \code{theta} is.
#' @seealso \code{\link{linkderiv}}, which routes to this generic by order.
#' @examples
#' d3linkfun(logit_link(), 0.5)
#' d3linkfun(log_link(), c(1, 2))
#' @export
d3linkfun <- S7::new_generic("d3linkfun", "x", fun = function(x, theta) S7::S7_dispatch())

#' @title 4th Derivative of Link Function
#' @description
#' The fourth derivative of \eqn{g(\theta)} with respect to \eqn{\theta}.
#' @param x An object of class \code{link}.
#' @param theta A numeric vector.
#' @return A numeric vector of the same length as \code{theta}, missing wherever
#'   \code{theta} is.
#' @seealso \code{\link{linkderiv}}, which routes to this generic by order.
#' @examples
#' d4linkfun(logit_link(), 0.5)
#' d4linkfun(log_link(), c(1, 2))
#' @export
d4linkfun <- S7::new_generic("d4linkfun", "x", fun = function(x, theta) S7::S7_dispatch())

#' @title 1st Derivative of Inverse Link Function
#' @description
#' The first derivative of \eqn{g^{-1}(\eta)} with respect to \eqn{\eta}.
#' @param x An object of class \code{link}.
#' @param eta A numeric vector.
#' @return A numeric vector of the same length as \code{eta}, missing wherever
#'   \code{eta} is.
#' @details
#' This and its higher-order siblings are the generics a modelling routine
#' working on the unconstrained scale actually wants. Call them directly rather
#' than through \code{\link{linkinvderiv}} in a hot loop: the router dispatches
#' once on itself and then again on the order-specific generic, which is about a
#' third of the cost of the call.
#' @seealso \code{\link{linkinvderiv}}, which routes to this generic by order.
#' @examples
#' dlinkinv(logit_link(), 0)      # p(1 - p) at p = 0.5
#' dlinkinv(log_link(), c(0, 1))
#' @export
dlinkinv <- S7::new_generic("dlinkinv", "x", fun = function(x, eta) S7::S7_dispatch())

#' @title 2nd Derivative of Inverse Link Function
#' @description
#' The second derivative of \eqn{g^{-1}(\eta)} with respect to \eqn{\eta}.
#' @param x An object of class \code{link}.
#' @param eta A numeric vector.
#' @return A numeric vector of the same length as \code{eta}, missing wherever
#'   \code{eta} is.
#' @seealso \code{\link{linkinvderiv}}, which routes to this generic by order.
#' @examples
#' d2linkinv(logit_link(), 0)     # zero, by symmetry about eta = 0
#' d2linkinv(probit_link(), 1)
#' @export
d2linkinv <- S7::new_generic("d2linkinv", "x", fun = function(x, eta) S7::S7_dispatch())

#' @title 3rd Derivative of Inverse Link Function
#' @description
#' The third derivative of \eqn{g^{-1}(\eta)} with respect to \eqn{\eta}.
#' @param x An object of class \code{link}.
#' @param eta A numeric vector.
#' @return A numeric vector of the same length as \code{eta}, missing wherever
#'   \code{eta} is.
#' @seealso \code{\link{linkinvderiv}}, which routes to this generic by order.
#' @examples
#' d3linkinv(logit_link(), 0)
#' d3linkinv(probit_link(), 1)
#' @export
d3linkinv <- S7::new_generic("d3linkinv", "x", fun = function(x, eta) S7::S7_dispatch())

#' @title 4th Derivative of Inverse Link Function
#' @description
#' The fourth derivative of \eqn{g^{-1}(\eta)} with respect to \eqn{\eta}.
#' @param x An object of class \code{link}.
#' @param eta A numeric vector.
#' @return A numeric vector of the same length as \code{eta}, missing wherever
#'   \code{eta} is.
#' @seealso \code{\link{linkinvderiv}}, which routes to this generic by order.
#' @examples
#' d4linkinv(logit_link(), 0)
#' d4linkinv(probit_link(), 1)
#' @export
d4linkinv <- S7::new_generic("d4linkinv", "x", fun = function(x, eta) S7::S7_dispatch())

#' @title Evaluate Derivative of Link Function by Order
#' @description
#' A convenience router over \code{\link{linkfun}} and the four
#' \code{d*linkfun} generics.
#' @param x An object of class \code{link}.
#' @param theta A numeric vector.
#' @param order An integer specifying the derivative order (0 to 4). Order 0 is
#'   the link function itself.
#' @return A numeric vector of the same length as \code{theta}: the requested
#'   derivative of \eqn{g} evaluated there.
#' @details
#' This dispatches twice, once on itself and once on the order-specific generic.
#' Where that matters, call \code{\link{dlinkfun}} and its siblings directly.
#' @examples
#' lk <- logit_link()
#' linkderiv(lk, 0.5, order = 0)   # the link itself
#' linkderiv(lk, 0.5, order = 1)
#'
#' # every order at once
#' vapply(0:4, function(k) linkderiv(lk, 0.25, order = k), numeric(1))
#' @export
linkderiv <- S7::new_generic("linkderiv", "x", fun = function(x, theta, order = 1) S7::S7_dispatch())

#' @title Evaluate Derivative of Inverse Link Function by Order
#' @description
#' A convenience router over \code{\link{linkinv}} and the four
#' \code{d*linkinv} generics.
#' @param x An object of class \code{link}.
#' @param eta A numeric vector.
#' @param order An integer specifying the derivative order (0 to 4). Order 0 is
#'   the inverse link itself.
#' @return A numeric vector of the same length as \code{eta}: the requested
#'   derivative of \eqn{g^{-1}} evaluated there.
#' @details
#' This dispatches twice, once on itself and once on the order-specific generic.
#' Where that matters, call \code{\link{dlinkinv}} and its siblings directly.
#' @examples
#' lk <- logit_link()
#' linkinvderiv(lk, 0, order = 0)  # the inverse link itself
#' linkinvderiv(lk, 0, order = 1)
#'
#' vapply(0:4, function(k) linkinvderiv(lk, 0.5, order = k), numeric(1))
#' @export
linkinvderiv <- S7::new_generic("linkinvderiv", "x", fun = function(x, eta, order = 1) S7::S7_dispatch())

#' @title Validate and Check a Link Object
#' @param x An object of class \code{link}.
#' @param tolerance Numeric tolerance for floating-point comparisons.
#' @param ... Additional arguments passed to methods.
#' @return Invisibly, a named list of check results; see
#'   \code{\link{check_link.link}} for its shape. Called mainly for the summary
#'   printed to the console.
#' @examples
#' check_link(sqrt_link())
#' @export
check_link <- S7::new_generic("check_link", "x", fun = function(x, tolerance = 1e-5, ...) S7::S7_dispatch())

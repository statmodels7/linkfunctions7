#' @title S7 Class for the Identity Link
#'
#' @description
#' Carries the identity transformation \eqn{\eta = \theta} on the whole real
#' line, for a parameter that needs no chart because it is already
#' unconstrained. Every derivative is a constant: the first is one and the
#' higher ones are zero, in both directions.
#'
#' [bounded_link()] returns an object of this class when neither endpoint is
#' given.
#'
#' @param link_name A character string naming the link, set by the
#'   constructor and shown by `print()`.
#' @param link_bounds A length-two numeric vector, the open interval in
#'   which the parameter lies. Set by the constructor; see Value for this link's.
#' @param link_params A list of the link's own parameters, empty where it has
#'   none. Set by the constructor.
#'
#' @return An S7 object of class `IdentityLink`, inheriting from [link()] and
#'   carrying its three properties `link_name`, `link_bounds` and
#'   `link_params`. Its `link_bounds` are `c(-Inf, Inf)` and it carries no link parameters.
#'
#' @section Methods:
#' Twelve methods are registered on this class: [linkfun()] and [linkinv()],
#' and the five derivative orders in each direction, [dlinkfun()] through
#' [d5linkfun()] going out and [dlinkinv()] through [d5linkinv()] coming
#' back. `linkfun()` and `linkinv()` return their argument. The first
#' derivative is one and every higher one is zero, each built by
#' [const_like()] so that the result is `NA` wherever the argument is `NA`,
#' although the derivative does not depend on it.
#'
#' @aliases linkfun.IdentityLink
#' @aliases linkinv.IdentityLink
#' @aliases dlinkfun.IdentityLink
#' @aliases d2linkfun.IdentityLink
#' @aliases d3linkfun.IdentityLink
#' @aliases d4linkfun.IdentityLink
#' @aliases d5linkfun.IdentityLink
#' @aliases dlinkinv.IdentityLink
#' @aliases d2linkinv.IdentityLink
#' @aliases d3linkinv.IdentityLink
#' @aliases d4linkinv.IdentityLink
#' @aliases d5linkinv.IdentityLink
#'
#' @seealso [identity_link()], the constructor users call.
#' @keywords internal
IdentityLink <- S7::new_class(
  name = "IdentityLink",
  parent = link
)

# --- Methods for IdentityLink ---

# Forward and inverse link functions
S7::method(linkfun, IdentityLink) <- function(x, theta) theta
S7::method(linkinv, IdentityLink) <- function(x, eta) eta

# Exact analytical derivatives of the link function (wrt theta)
S7::method(dlinkfun, IdentityLink) <- function(x, theta) const_like(theta, 1)
S7::method(d2linkfun, IdentityLink) <- function(x, theta) const_like(theta, 0)
S7::method(d3linkfun, IdentityLink) <- function(x, theta) const_like(theta, 0)
S7::method(d4linkfun, IdentityLink) <- function(x, theta) const_like(theta, 0)
S7::method(d5linkfun, IdentityLink) <- function(x, theta) const_like(theta, 0)

# Exact analytical derivatives of the inverse link function (wrt eta)
S7::method(dlinkinv, IdentityLink) <- function(x, eta) const_like(eta, 1)
S7::method(d2linkinv, IdentityLink) <- function(x, eta) const_like(eta, 0)
S7::method(d3linkinv, IdentityLink) <- function(x, eta) const_like(eta, 0)
S7::method(d4linkinv, IdentityLink) <- function(x, eta) const_like(eta, 0)
S7::method(d5linkinv, IdentityLink) <- function(x, eta) const_like(eta, 0)

#' @title The Identity Link Function
#'
#' @include generics.R
#' @include link_class.R
#' @description
#' The identity link \eqn{\eta = \theta}, for a parameter that is already
#' unconstrained.
#' @details
#' The identity link is \eqn{\eta = \theta}, and its inverse is
#' \eqn{\theta = \eta}.
#'
#' The first derivative is 1 in both directions, and every higher derivative
#' is 0.
#'
#' The domain of \eqn{\theta} is the whole real line.
#'
#' @return An S7 object of class `IdentityLink`, inheriting from [link()], whose
#'   methods compute the link, its inverse and their derivatives to the fifth
#'   order.
#'
#' @examples
#' lk <- identity_link()
#' lk
#'
#' linkfun(lk, c(-1, 0, 1))
#' linkinv(lk, c(-1, 0, 1))
#'
#' # the first derivative is 1 and every higher one is 0 ...
#' dlinkfun(lk, c(-1, 0, 1))
#' d2linkfun(lk, c(-1, 0, 1))
#'
#' # ... but missingness is still propagated, not swallowed by the constant
#' dlinkfun(lk, c(1, NA))
#'
#' @seealso [link()]
#' @export
identity_link <- function() {
  IdentityLink(
    link_name = "identity",
    link_bounds = c(-Inf, Inf),
    
    # The identity link requires no additional mathematical parameters
    link_params = NULL
  )
}

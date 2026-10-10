#' @title S7 Class for the Sqrt Link
#'
#' @description
#' Carries the square-root transformation \eqn{\eta = \sqrt{\theta}} on
#' \eqn{(0, \infty)}, with inverse \eqn{\theta = \eta^2}.
#'
#' Its image is \eqn{(0, \infty)}. On the whole real line the inverse link is
#' not one-to-one, since \eqn{\eta} and \eqn{-\eta} give the same
#' \eqn{\theta}; [eta_bounds()] returns the range on which it is.
#'
#' @param link_name A character string naming the link, set by the
#'   constructor and shown by `print()`.
#' @param link_bounds A length-two numeric vector, the open interval in
#'   which the parameter lies. Set by the constructor; see Value for this link's.
#' @param link_params A list of the link's own parameters, empty where it has
#'   none. Set by the constructor.
#'
#' @return An S7 object of class `SqrtLink`, inheriting from [link()] and
#'   carrying its three properties `link_name`, `link_bounds` and
#'   `link_params`. Its `link_bounds` are `c(0, Inf)` and it carries no link parameters.
#'
#' @section Methods:
#' Twelve methods are registered on this class: [linkfun()] and [linkinv()],
#' and the five derivative orders in each direction, [dlinkfun()] through
#' [d5linkfun()] going out and [dlinkinv()] through [d5linkinv()] coming
#' back. The forward derivatives are half-integer falling factorials. The
#' inverse map is \eqn{\eta^2}, so its derivatives terminate: the second is
#' the constant two and the third to fifth are exactly zero, all built by
#' [const_like()] so that a missing value still propagates.
#'
#' @aliases linkfun.SqrtLink
#' @aliases linkinv.SqrtLink
#' @aliases dlinkfun.SqrtLink
#' @aliases d2linkfun.SqrtLink
#' @aliases d3linkfun.SqrtLink
#' @aliases d4linkfun.SqrtLink
#' @aliases d5linkfun.SqrtLink
#' @aliases dlinkinv.SqrtLink
#' @aliases d2linkinv.SqrtLink
#' @aliases d3linkinv.SqrtLink
#' @aliases d4linkinv.SqrtLink
#' @aliases d5linkinv.SqrtLink
#'
#' @seealso [sqrt_link()], the constructor users call.
#' @keywords internal
SqrtLink <- S7::new_class(
  name = "SqrtLink",
  parent = link
)

# --- Methods for SqrtLink ---

# Forward and inverse link functions
S7::method(linkfun, SqrtLink) <- function(x, theta) sqrt(theta)
S7::method(linkinv, SqrtLink) <- function(x, eta) eta^2

# Exact analytical derivatives of the link function (wrt theta)
S7::method(dlinkfun, SqrtLink) <- function(x, theta) {
  1 / (2 * sqrt(theta))
}
S7::method(d2linkfun, SqrtLink) <- function(x, theta) {
  -1 / (4 * (theta^1.5))
}
S7::method(d3linkfun, SqrtLink) <- function(x, theta) {
  3 / (8 * (theta^2.5))
}
S7::method(d4linkfun, SqrtLink) <- function(x, theta) {
  -15 / (16 * (theta^3.5))
}
S7::method(d5linkfun, SqrtLink) <- function(x, theta) {
  105 / (32 * (theta^4.5))
}

# Exact analytical derivatives of the inverse link function (wrt eta)
# 3rd and 4th derivatives uniquely vanish to exactly 0 for this quadratic form.
S7::method(dlinkinv, SqrtLink) <- function(x, eta) 2 * eta
S7::method(d2linkinv, SqrtLink) <- function(x, eta) const_like(eta, 2)
S7::method(d3linkinv, SqrtLink) <- function(x, eta) const_like(eta, 0)
S7::method(d4linkinv, SqrtLink) <- function(x, eta) const_like(eta, 0)
S7::method(d5linkinv, SqrtLink) <- function(x, eta) const_like(eta, 0)

#' @title The Square Root Link Function
#'
#' @include generics.R
#' @include link_class.R
#' @description
#' The square-root link \eqn{\eta = \sqrt{\theta}} on \eqn{(0, \infty)};
#' its image is \eqn{(0, \infty)}.
#' @details
#' The square-root link is \eqn{\eta = \sqrt{\theta}}, with inverse
#' \eqn{\theta = \eta^2}.
#'
#' The inverse \eqn{\eta^2} is defined for a negative \eqn{\eta} as well, but
#' there it is not one-to-one, since \eqn{\eta} and \eqn{-\eta} give the same
#' \eqn{\theta}. The linear predictor is therefore meant to stay positive;
#' [eta_bounds()] returns this range.
#'
#' The domain of \eqn{\theta} is \eqn{(0, \infty)}.
#'
#' @return An S7 object of class `SqrtLink`, inheriting from [link()], whose
#'   methods compute the link, its inverse and their derivatives to the fifth
#'   order.
#'
#' @examples
#' lk <- sqrt_link()
#' lk
#'
#' theta <- c(0.25, 1, 4)
#' eta <- linkfun(lk, theta)
#' eta
#' linkinv(lk, eta)
#'
#' # the inverse is a quadratic, so the third and fourth derivatives vanish
#' d2linkinv(lk, c(1, 2))
#' d3linkinv(lk, c(1, 2))
#'
#' # the same link as the power family at lambda = 1/2
#' linkfun(power_link(0.5), 4)
#'
#' @seealso [link()], [power_link()], [log_link()]
#' @export
sqrt_link <- function() {
  SqrtLink(
    link_name = "sqrt",
    link_bounds = c(0, Inf),
    
    # The standard square root link requires no additional parameters
    link_params = NULL
  )
}

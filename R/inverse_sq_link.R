#' @title S7 Class for the InverseSq Link
#'
#' @description
#' Carries the inverse-square transformation \eqn{\eta = 1/\theta^2} on
#' \eqn{(0, \infty)}, with inverse \eqn{\theta = 1/\sqrt{\eta}}.
#'
#' It is the canonical link of the inverse Gaussian family. Its image is
#' \eqn{(0, \infty)}, so the inverse link is defined only for a positive linear
#' predictor, and like [inverse_link()] the map is decreasing.
#'
#' @param link_name A character string naming the link, set by the
#'   constructor and shown by `print()`.
#' @param link_bounds A length-two numeric vector, the open interval in
#'   which the parameter lies. Set by the constructor; see Value for this link's.
#' @param link_params A list of the link's own parameters, empty where it has
#'   none. Set by the constructor.
#'
#' @return An S7 object of class `InverseSqLink`, inheriting from [link()] and
#'   carrying its three properties `link_name`, `link_bounds` and
#'   `link_params`. Its `link_bounds` are `c(0, Inf)` and it carries no link parameters.
#'
#' @section Methods:
#' Twelve methods are registered on this class: [linkfun()] and [linkinv()],
#' and the five derivative orders in each direction, [dlinkfun()] through
#' [d5linkfun()] going out and [dlinkinv()] through [d5linkinv()] coming
#' back. The two directions are not the same map, \eqn{1/\theta^2} going
#' forward and \eqn{\eta^{-1/2}} coming back, so the two sets of derivatives
#' are written out separately, with integer falling factorials in the forward
#' direction and half-integer powers in the inverse one.
#'
#' @aliases linkfun.InverseSqLink
#' @aliases linkinv.InverseSqLink
#' @aliases dlinkfun.InverseSqLink
#' @aliases d2linkfun.InverseSqLink
#' @aliases d3linkfun.InverseSqLink
#' @aliases d4linkfun.InverseSqLink
#' @aliases d5linkfun.InverseSqLink
#' @aliases dlinkinv.InverseSqLink
#' @aliases d2linkinv.InverseSqLink
#' @aliases d3linkinv.InverseSqLink
#' @aliases d4linkinv.InverseSqLink
#' @aliases d5linkinv.InverseSqLink
#'
#' @seealso [inverse_sq_link()], the constructor users call.
#' @keywords internal
InverseSqLink <- S7::new_class(
  name = "InverseSqLink",
  parent = link
)

# --- Methods for InverseSqLink ---

# Forward and inverse link functions
S7::method(linkfun, InverseSqLink) <- function(x, theta) 1 / (theta^2)
S7::method(linkinv, InverseSqLink) <- function(x, eta) 1 / sqrt(eta)

# Exact analytical derivatives of the link function (wrt theta)
S7::method(dlinkfun, InverseSqLink) <- function(x, theta) -2 / (theta^3)
S7::method(d2linkfun, InverseSqLink) <- function(x, theta)  6 / (theta^4)
S7::method(d3linkfun, InverseSqLink) <- function(x, theta) -24 / (theta^5)
S7::method(d4linkfun, InverseSqLink) <- function(x, theta) 120 / (theta^6)
S7::method(d5linkfun, InverseSqLink) <- function(x, theta) -720 / (theta^7)

# Exact analytical derivatives of the inverse link function (wrt eta)
# Utilizing explicit fractions to optimize numeric evaluations
S7::method(dlinkinv, InverseSqLink) <- function(x, eta) -1 / (2 * eta^1.5)
S7::method(d2linkinv, InverseSqLink) <- function(x, eta)  3 / (4 * eta^2.5)
S7::method(d3linkinv, InverseSqLink) <- function(x, eta) -15 / (8 * eta^3.5)
S7::method(d4linkinv, InverseSqLink) <- function(x, eta) 105 / (16 * eta^4.5)
S7::method(d5linkinv, InverseSqLink) <- function(x, eta) -945 / (32 * eta^5.5)

#' @title The Inverse Square Link Function
#'
#' @include generics.R
#' @include link_class.R
#' @description
#' The inverse-square link \eqn{\eta = 1/\theta^2} on \eqn{(0, \infty)},
#' the canonical link of the inverse Gaussian family; its image is
#' \eqn{(0, \infty)}.
#' @details
#' The inverse-square link is \eqn{\eta = 1/\theta^2}, with inverse
#' \eqn{\theta = 1/\sqrt{\eta}}.
#'
#' It is the canonical link of the inverse Gaussian family, whose variance is
#' proportional to the cube of the mean.
#'
#' The domain of \eqn{\theta} is \eqn{(0, \infty)}, and so is the image of
#' the link. The inverse link and its derivatives return `NaN` for a negative
#' \eqn{\eta}, so the linear predictor must stay positive during
#' optimization; [eta_bounds()] returns this range.
#'
#' @return An S7 object of class `InverseSqLink`, inheriting from [link()], whose
#'   methods compute the link, its inverse and their derivatives to the fifth
#'   order.
#'
#' @examples
#' lk <- inverse_sq_link()
#' lk
#'
#' theta <- c(0.5, 1, 2)
#' eta <- linkfun(lk, theta)  # 1 / theta^2
#' eta
#' linkinv(lk, eta)
#'
#' # the canonical link of the inverse Gaussian; eta must stay positive
#' dlinkinv(lk, c(0.5, 1, 4))
#'
#' @seealso [link()], [inverse_link()]
#' @export
inverse_sq_link <- function() {
  InverseSqLink(
    link_name = "inverse_sq",
    link_bounds = c(0, Inf),
    
    # The inverse square link requires no additional mathematical parameters
    link_params = NULL
  )
}

#' @title S7 Class for the Inverse Link
#'
#' @description
#' Carries the reciprocal transformation \eqn{\eta = 1/\theta} on
#' \eqn{(0, \infty)}, which is its own inverse.
#'
#' It is the canonical link of the Gamma family. Its image is
#' \eqn{(0, \infty)}, so the inverse link is defined only for a positive linear
#' predictor, and the map is decreasing: a larger parameter gives a smaller
#' \eqn{\eta}.
#'
#' @param link_name A character string naming the link, set by the
#'   constructor and shown by `print()`.
#' @param link_bounds A length-two numeric vector, the open interval in
#'   which the parameter lies. Set by the constructor; see Value for this link's.
#' @param link_params A list of the link's own parameters, empty where it has
#'   none. Set by the constructor.
#'
#' @return An S7 object of class `InverseLink`, inheriting from [link()] and
#'   carrying its three properties `link_name`, `link_bounds` and
#'   `link_params`. Its `link_bounds` are `c(0, Inf)` and it carries no link parameters.
#'
#' @section Methods:
#' Twelve methods are registered on this class: [linkfun()] and [linkinv()],
#' and the five derivative orders in each direction, [dlinkfun()] through
#' [d5linkfun()] going out and [dlinkinv()] through [d5linkinv()] coming
#' back. The link is its own inverse, so the two directions carry the same
#' expressions: \eqn{1/\theta} and \eqn{1/\eta}, with the five derivatives
#' \eqn{(-1)^k k!\,z^{-(k+1)}} written out in both.
#'
#' @aliases linkfun.InverseLink
#' @aliases linkinv.InverseLink
#' @aliases dlinkfun.InverseLink
#' @aliases d2linkfun.InverseLink
#' @aliases d3linkfun.InverseLink
#' @aliases d4linkfun.InverseLink
#' @aliases d5linkfun.InverseLink
#' @aliases dlinkinv.InverseLink
#' @aliases d2linkinv.InverseLink
#' @aliases d3linkinv.InverseLink
#' @aliases d4linkinv.InverseLink
#' @aliases d5linkinv.InverseLink
#'
#' @seealso [inverse_link()], the constructor users call.
#' @keywords internal
InverseLink <- S7::new_class(
  name = "InverseLink",
  parent = link
)

# --- Methods for InverseLink ---

# Forward and inverse link functions
S7::method(linkfun, InverseLink) <- function(x, theta) 1 / theta
S7::method(linkinv, InverseLink) <- function(x, eta) 1 / eta

# Exact analytical derivatives of the link function (wrt theta)
S7::method(dlinkfun, InverseLink) <- function(x, theta) -1 / (theta^2)
S7::method(d2linkfun, InverseLink) <- function(x, theta)  2 / (theta^3)
S7::method(d3linkfun, InverseLink) <- function(x, theta) -6 / (theta^4)
S7::method(d4linkfun, InverseLink) <- function(x, theta) 24 / (theta^5)
S7::method(d5linkfun, InverseLink) <- function(x, theta) -120 / (theta^6)

# Exact analytical derivatives of the inverse link function (wrt eta)
# Due to the symmetric nature of f(x) = 1/x, these are structurally identical
# to the link function derivatives, but evaluated at eta.
S7::method(dlinkinv, InverseLink) <- function(x, eta) -1 / (eta^2)
S7::method(d2linkinv, InverseLink) <- function(x, eta)  2 / (eta^3)
S7::method(d3linkinv, InverseLink) <- function(x, eta) -6 / (eta^4)
S7::method(d4linkinv, InverseLink) <- function(x, eta) 24 / (eta^5)
S7::method(d5linkinv, InverseLink) <- function(x, eta) -120 / (eta^6)

#' @title The Inverse (Reciprocal) Link Function
#'
#' @include generics.R
#' @include link_class.R
#' @description
#' The reciprocal link \eqn{\eta = 1/\theta} on \eqn{(0, \infty)}, the
#' canonical link of the Gamma family; its image is \eqn{(0, \infty)},
#' not the whole real line.
#' @details
#' The inverse link is \eqn{\eta = 1/\theta}, and it is its own inverse,
#' \eqn{\theta = 1/\eta}.
#'
#' It is the canonical link of the Gamma family, used when the mean is
#' modeled as the reciprocal of a linear predictor.
#'
#' The domain of \eqn{\theta} is \eqn{(0, \infty)}, and so is the image of
#' the link, so the linear predictor must stay positive during optimization;
#' [eta_bounds()] returns this range.
#'
#' @return An S7 object of class `InverseLink`, inheriting from [link()], whose
#'   methods compute the link, its inverse and their derivatives to the fifth
#'   order.
#'
#' @examples
#' lk <- inverse_link()
#' lk
#'
#' theta <- c(0.5, 1, 2)
#' eta <- linkfun(lk, theta)
#' eta
#' linkinv(lk, eta)           # the map is its own inverse
#'
#' dlinkfun(lk, theta)
#'
#' # the canonical link for a Gamma mean; note eta must keep one sign
#' linkinv(lk, c(0.5, 2))
#'
#' @seealso [link()], [identity_link()]
#' @export
inverse_link <- function() {
  InverseLink(
    link_name = "inverse",
    link_bounds = c(0, Inf),
    
    # The inverse link requires no additional mathematical parameters
    link_params = NULL
  )
}

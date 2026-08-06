#' @title S7 Class for the Rhobit Link
#'
#' @description
#' The class \code{\link{rhobit_link}} instantiates.
#'
#' @return An S7 object of class \code{RhobitLink}, inheriting from
#'   \code{\link{link}}.
#'
#' @seealso \code{\link{rhobit_link}}, the constructor users call.
#' @keywords internal
RhobitLink <- S7::new_class(
  name = "RhobitLink",
  parent = link
)

# --- Methods for RhobitLink ---

# Forward and inverse link functions
S7::method(linkfun, RhobitLink) <- function(x, theta) atanh(theta)
S7::method(linkinv, RhobitLink) <- function(x, eta) tanh(eta)

# Exact analytical derivatives of the link function (wrt theta)
S7::method(dlinkfun, RhobitLink) <- function(x, theta) lk_rhobit_fwd_cpp(theta, 1L)
S7::method(d2linkfun, RhobitLink) <- function(x, theta) lk_rhobit_fwd_cpp(theta, 2L)
S7::method(d3linkfun, RhobitLink) <- function(x, theta) lk_rhobit_fwd_cpp(theta, 3L)
S7::method(d4linkfun, RhobitLink) <- function(x, theta) lk_rhobit_fwd_cpp(theta, 4L)

# Exact analytical derivatives of the inverse link function (wrt eta)
# All of these are evaluated as polynomials in t = tanh(eta)
S7::method(dlinkinv, RhobitLink) <- function(x, eta) lk_rhobit_inv_cpp(eta, 1L)
S7::method(d2linkinv, RhobitLink) <- function(x, eta) lk_rhobit_inv_cpp(eta, 2L)
S7::method(d3linkinv, RhobitLink) <- function(x, eta) lk_rhobit_inv_cpp(eta, 3L)
S7::method(d4linkinv, RhobitLink) <- function(x, eta) lk_rhobit_inv_cpp(eta, 4L)

#' @title The Rhobit (Fisher's z) Link Function
#'
#' @include generics.R
#' @include link_class.R
#' @description
#' The rhobit link
#' \eqn{\eta = \mathrm{atanh}(\theta) = \log((1+\theta)/(1-\theta))/2}
#' on \eqn{(-1, 1)}, Fisher's z; the natural link for a correlation.
#' @details
#' The Rhobit link is defined mathematically using the inverse hyperbolic tangent function:
#' \eqn{\eta = \text{arctanh}(\theta) = \frac{1}{2} \log\left(\frac{1 + \theta}{1 - \theta}\right)}.
#'
#' The inverse link is the hyperbolic tangent function:
#' \eqn{\theta = \tanh(\eta) = \frac{\exp(2\eta) - 1}{\exp(2\eta) + 1}}.
#'
#' The valid mathematical domain of \eqn{\theta} is exactly \code{c(-1, 1)}.
#'
#' @return An S7 object of class \code{RhobitLink} (inheriting from \code{link}) containing the transformation functions
#' and their exact analytical derivatives up to the fourth order.
#'
#' @examples
#' lk <- rhobit_link()
#' lk
#'
#' # built for parameters constrained to (-1, 1), such as a correlation
#' rho <- c(-0.9, 0, 0.9)
#' eta <- linkfun(lk, rho)    # Fisher's z
#' eta
#' linkinv(lk, eta)
#'
#' # the inverse is tanh, so its first derivative is the squared sech
#' dlinkinv(lk, 0)
#'
#' check_link(lk)
#'
#' @seealso \code{\link{link}}, \code{\link{logit_link}}
#' @export
rhobit_link <- function() {
  RhobitLink(
    link_name = "rhobit",
    link_bounds = c(-1, 1),
    link_params = NULL
  )
}

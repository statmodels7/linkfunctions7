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

# Exact analytical derivatives of the inverse link function (wrt eta)
# Utilizing explicit fractions to optimize numeric evaluations
S7::method(dlinkinv, InverseSqLink) <- function(x, eta) -1 / (2 * eta^1.5)
S7::method(d2linkinv, InverseSqLink) <- function(x, eta)  3 / (4 * eta^2.5)
S7::method(d3linkinv, InverseSqLink) <- function(x, eta) -15 / (8 * eta^3.5)
S7::method(d4linkinv, InverseSqLink) <- function(x, eta) 105 / (16 * eta^4.5)

#' @title The Inverse Square Link Function
#'
#' @include generics.R
#' @include link_class.R
#' @description
#' Creates an S7 object of class \code{link} implementing the Inverse Square transformation.
#' This is the canonical link function for the Inverse Gaussian distribution.
#'
#' @details
#' The Inverse Square link is defined mathematically as \eqn{\eta = 1 / \theta^2}.
#' Consequently, the inverse link function is derived as \eqn{\theta = 1 / \sqrt{\eta}}.
#'
#' This specific link function is predominantly utilized in Generalized Linear Models (GLMs)
#' assuming an Inverse Gaussian response distribution. In such frameworks, the variance
#' is proportional to the cube of the mean (\eqn{\text{Var}(Y) \propto \mu^3}).
#'
#' \strong{Domain and Optimization Constraints:}
#' Both the parameter \eqn{\theta} and the linear predictor \eqn{\eta} must be strictly
#' positive. The valid mathematical domain for \eqn{\theta} is \code{c(0, Inf)}. During
#' optimization routines (e.g., Fisher Scoring or Newton-Raphson), extreme care must be
#' taken to ensure the linear predictor \eqn{\eta > 0}. Evaluating the inverse link or
#' its derivatives at non-positive values of \eqn{\eta} will inevitably result in \code{NaN}s
#' due to fractional powers and square root operations.
#'
#' @return An S7 object of class \code{InverseSqLink} (inheriting from \code{link}) containing the transformation functions
#' and their exact analytical derivatives up to the fourth order.
#'
#' @seealso \code{\link{link}}, \code{\link{inverse_link}}
#' @export
inverse_sq_link <- function() {
  InverseSqLink(
    link_name = "inverse_sq",
    link_bounds = c(0, Inf),
    
    # The inverse square link requires no additional mathematical parameters
    link_params = NULL
  )
}
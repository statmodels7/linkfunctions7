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

# Exact analytical derivatives of the inverse link function (wrt eta)
S7::method(dlinkinv, IdentityLink) <- function(x, eta) const_like(eta, 1)
S7::method(d2linkinv, IdentityLink) <- function(x, eta) const_like(eta, 0)
S7::method(d3linkinv, IdentityLink) <- function(x, eta) const_like(eta, 0)
S7::method(d4linkinv, IdentityLink) <- function(x, eta) const_like(eta, 0)

#' @title The Identity Link Function
#'
#' @include generics.R
#' @include link_class.R
#' @description
#' Creates an S7 object of class \code{link} implementing the Identity transformation.
#' This is the canonical link function for the mean parameter of the Normal (Gaussian) 
#' distribution.
#'
#' @details
#' The Identity link is defined simply as \eqn{\eta = \theta}.
#' Consequently, the inverse link is also \eqn{\theta = \eta}.
#'
#' All first derivatives are constant (equal to 1), and all higher-order derivatives
#' up to the fourth order are exactly zero.
#'
#' The domain of \eqn{\theta} is unbounded, meaning the valid domain is \code{c(-Inf, Inf)}.
#'
#' @return An S7 object of class \code{IdentityLink} (inheriting from \code{link}) containing the transformation functions
#' and their exact analytical derivatives up to the fourth order.
#'
#' @seealso \code{\link{link}}
#' @export
identity_link <- function() {
  IdentityLink(
    link_name = "identity",
    link_bounds = c(-Inf, Inf),
    
    # The identity link requires no additional mathematical parameters
    link_params = NULL
  )
}
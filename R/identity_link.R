#' @title The Identity Link Function
#'
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
#' @return An S7 object of class \code{link} containing the transformation functions
#' and their exact analytical derivatives up to the fourth order.
#'
#' @seealso \code{\link{link}}
#'
#' @export
identity_link <- function() {
  link(
    link_name = "identity",
    link_bounds = c(-Inf, Inf),
    
    # The identity link requires no additional mathematical parameters
    link_params = NULL,
    
    # Forward and inverse link functions
    linkfun = function(theta) theta,
    linkinv = function(eta) eta,
    
    # Exact analytical derivatives of the link function (wrt theta)
    dlinkfun  = function(theta) rep(1, length(theta)),
    d2linkfun = function(theta) rep(0, length(theta)),
    d3linkfun = function(theta) rep(0, length(theta)),
    d4linkfun = function(theta) rep(0, length(theta)),
    
    # Exact analytical derivatives of the inverse link function (wrt eta)
    dlinkinv  = function(eta) rep(1, length(eta)),
    d2linkinv = function(eta) rep(0, length(eta)),
    d3linkinv = function(eta) rep(0, length(eta)),
    d4linkinv = function(eta) rep(0, length(eta))
  )
}
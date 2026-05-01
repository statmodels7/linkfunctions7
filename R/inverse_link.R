#' @title The Inverse (Reciprocal) Link Function
#'
#' @description
#' Creates an S7 object of class \code{link} implementing the reciprocal transformation.
#'
#' @details
#' The Inverse link is defined as \eqn{\eta = 1/\theta}.
#' The inverse link function is therefore perfectly symmetric: \eqn{\theta = 1/\eta}.
#'
#' This link is typically used for modeling positive continuous data where the mean is 
#' inversely proportional to the linear predictor (e.g., in Gamma regression).
#'
#' The domain of \eqn{\theta} is conventionally \code{c(0, Inf)}. Care must be taken 
#' to ensure the linear predictor \eqn{\eta} remains strictly positive (or strictly 
#' negative) during optimization to avoid division by zero or mapping to invalid 
#' negative parameter values.
#'
#' @return An S7 object of class \code{link} containing the transformation functions
#' and their exact analytical derivatives up to the fourth order.
#'
#' @seealso \code{\link{link}}, \code{\link{identity_link}}
#'
#' @export
inverse_link <- function() {
  link(
    link_name = "inverse",
    link_bounds = c(0, Inf),
    
    # The inverse link requires no additional mathematical parameters
    link_params = NULL,
    
    # Forward and inverse link functions
    linkfun = function(theta) 1 / theta,
    linkinv = function(eta) 1 / eta,
    
    # Exact analytical derivatives of the link function (wrt theta)
    dlinkfun  = function(theta) -1 / (theta^2),
    d2linkfun = function(theta)  2 / (theta^3),
    d3linkfun = function(theta) -6 / (theta^4),
    d4linkfun = function(theta) 24 / (theta^5),
    
    # Exact analytical derivatives of the inverse link function (wrt eta)
    # Due to the symmetric nature of f(x) = 1/x, these are structurally identical
    # to the link function derivatives, but evaluated at eta.
    dlinkinv  = function(eta) -1 / (eta^2),
    d2linkinv = function(eta)  2 / (eta^3),
    d3linkinv = function(eta) -6 / (eta^4),
    d4linkinv = function(eta) 24 / (eta^5)
  )
}
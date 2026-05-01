#' @title The Power Link Function
#'
#' @description
#' Creates an S7 object of class \code{link} implementing the Power transformation family.
#' This function generates a specific link function based on the provided power parameter \code{lambda}.
#' It elegantly includes the special case where \code{lambda = 0}, which dynamically returns the Log link.
#'
#' @param lambda A numeric value defining the power of the transformation. Defaults to 1.
#'
#' @details
#' The Power link is defined mathematically as \eqn{\eta = \theta^\lambda}.
#' Consequently, the inverse link is derived as \eqn{\theta = \eta^{1/\lambda}}.
#'
#' \strong{Special Case (Box-Cox continuity):}
#' If \code{lambda = 0}, the function mathematically approaches \eqn{\log(\theta)}. 
#' In this scenario, the function automatically instantiates and returns a \code{\link{log_link}} 
#' object, modifying its internal state to reflect the \code{lambda = 0} parameter.
#'
#' Common special cases include:
#' \itemize{
#'   \item \code{lambda = 1}: Identity link.
#'   \item \code{lambda = 0.5}: Square-root link.
#'   \item \code{lambda = -1}: Inverse link.
#'   \item \code{lambda = 0}: Log link.
#' }
#'
#' The mathematical domain of \eqn{\theta} is \code{c(0, Inf)}. Depending on the value 
#' of \code{lambda}, extreme care must be taken during numerical optimization to guarantee 
#' that \eqn{\eta} remains strictly positive to avoid \code{NaN}s from fractional exponents.
#'
#' @return An S7 object of class \code{link} containing the transformation functions, 
#' their exact analytical derivatives up to the fourth order, and the \code{link_params} 
#' property storing the value of \code{lambda}.
#'
#' @seealso \code{\link{link}}, \code{\link{log_link}}, \code{\link{identity_link}}
#'
#' @export
power_link <- function(lambda = 1) {
  # Handle Box-Cox continuity limit utilizing the existing log_link
  if (lambda == 0) {
    o <- log_link()
    o@link_name <- "power(lambda=0)"
    o@link_params <- list(lambda = 0)
    return(o)
  }

  k <- 1 / lambda

  link(
    link_name = paste0("power(lambda=", round(lambda, 5), ")"),
    link_bounds = c(0, Inf),
    link_params = list(lambda = lambda),
    
    linkfun = function(theta) theta^lambda,
    linkinv = function(eta) eta^k,
    
    dlinkfun  = function(theta) lambda * (theta^(lambda - 1)),
    d2linkfun = function(theta) lambda * (lambda - 1) * (theta^(lambda - 2)),
    d3linkfun = function(theta) lambda * (lambda - 1) * (lambda - 2) * (theta^(lambda - 3)),
    d4linkfun = function(theta) lambda * (lambda - 1) * (lambda - 2) * (lambda - 3) * (theta^(lambda - 4)),
    
    dlinkinv  = function(eta) k * (eta^(k - 1)),
    d2linkinv = function(eta) k * (k - 1) * (eta^(k - 2)),
    d3linkinv = function(eta) k * (k - 1) * (k - 2) * (eta^(k - 3)),
    d4linkinv = function(eta) k * (k - 1) * (k - 2) * (k - 3) * (eta^(k - 4))
  )
}
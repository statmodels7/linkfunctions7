#' @title S7 Class for Statistical Link Functions
#'
#' @import S7
#' @description
#' A strictly typed S7 object that encapsulates all mathematical properties
#' of a link function. It contains the transformation function, the inverse
#' transformation, and the exact analytical derivatives up to the fourth order
#' for both directions.
#'
#' @details
#' Objects of class \code{link} are instantiated using the S7 object system.
#'
#' The object assumes the following mathematical notation:
#' \itemize{
#'   \item \eqn{\theta}: The response parameter (e.g., probability, mean, dispersion).
#'   \item \eqn{\eta}: The linear predictor (unconstrained scale).
#' }
#' The relationship is defined as \eqn{\eta = g(\theta)} (link function) and
#' \eqn{\theta = g^{-1}(\eta)} (inverse link function).
#'
#' @param link_name A character string identifying the link (e.g., "logit").
#' @param link_bounds A numeric vector of length 2 \code{c(lower, upper)} defining the valid domain for \eqn{\theta}.
#' @param link_params A list or vector of additional parameters required to define the link, or \code{NULL}.
#' @param linkfun Function \eqn{g(\theta) \rightarrow \eta}.
#' @param linkinv Function \eqn{g^{-1}(\eta) \rightarrow \theta}.
#' @param dlinkfun First derivative of the link function.
#' @param d2linkfun Second derivative of the link function.
#' @param d3linkfun Third derivative of the link function.
#' @param d4linkfun Fourth derivative of the link function.
#' @param dlinkinv First derivative of the inverse link function.
#' @param d2linkinv Second derivative of the inverse link function.
#' @param d3linkinv Third derivative of the inverse link function.
#' @param d4linkinv Fourth derivative of the inverse link function.
#'
#' @export
link <- S7::new_class(
  name = "link",
  properties = list(
    link_name = S7::class_character,
    link_bounds = S7::class_numeric,
    link_params = S7::class_any,
    linkfun = S7::class_function,
    linkinv = S7::class_function,
    
    # Derivatives of the link function with respect to theta
    dlinkfun = S7::class_function,
    d2linkfun = S7::class_function,
    d3linkfun = S7::class_function,
    d4linkfun = S7::class_function,
    
    # Derivatives of the inverse link function with respect to eta
    dlinkinv = S7::class_function,
    d2linkinv = S7::class_function,
    d3linkinv = S7::class_function,
    d4linkinv = S7::class_function
  ),
  
  validator = function(self) {
    # Ensure bounds contain exactly two numeric elements
    if (length(self@link_bounds) != 2) {
      return("Property 'link_bounds' must be a numeric vector of length 2: c(lower, upper).")
    }
    
    # Ensure logical domain definition
    if (self@link_bounds[1] >= self@link_bounds[2]) {
      return("The lower bound must be strictly less than the upper bound.")
    }
  }
)

PowerLink <- S7::new_class(
  name = "PowerLink",
  parent = link,
  properties = list(
    lambda = S7::class_numeric
  )
)

# --- Methods for PowerLink ---

S7::method(linkfun, PowerLink) <- function(x, theta) theta^x@lambda
S7::method(linkinv, PowerLink) <- function(x, eta) eta^(1 / x@lambda)

S7::method(dlinkfun, PowerLink) <- function(x, theta) x@lambda * (theta^(x@lambda - 1))
S7::method(d2linkfun, PowerLink) <- function(x, theta) x@lambda * (x@lambda - 1) * (theta^(x@lambda - 2))
S7::method(d3linkfun, PowerLink) <- function(x, theta) x@lambda * (x@lambda - 1) * (x@lambda - 2) * (theta^(x@lambda - 3))
S7::method(d4linkfun, PowerLink) <- function(x, theta) x@lambda * (x@lambda - 1) * (x@lambda - 2) * (x@lambda - 3) * (theta^(x@lambda - 4))

S7::method(dlinkinv, PowerLink) <- function(x, eta) { k <- 1 / x@lambda; k * (eta^(k - 1)) }
S7::method(d2linkinv, PowerLink) <- function(x, eta) { k <- 1 / x@lambda; k * (k - 1) * (eta^(k - 2)) }
S7::method(d3linkinv, PowerLink) <- function(x, eta) { k <- 1 / x@lambda; k * (k - 1) * (k - 2) * (eta^(k - 3)) }
S7::method(d4linkinv, PowerLink) <- function(x, eta) { k <- 1 / x@lambda; k * (k - 1) * (k - 2) * (k - 3) * (eta^(k - 4)) }

#' @title The Power Link Function
#'
#' @include generics.R
#' @include link_class.R
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
#' @return An S7 object of class \code{PowerLink} (inheriting from \code{link}), 
#' or an object of class \code{LogLink} if \code{lambda = 0}.
#'
#' @seealso \code{\link{link}}, \code{\link{log_link}}, \code{\link{identity_link}}
#' @export
power_link <- function(lambda = 1) {
  # Handle Box-Cox continuity limit utilizing the existing log_link
  if (lambda == 0) {
    o <- log_link()
    o@link_name <- "power(lambda=0)"
    o@link_params <- list(lambda = 0)
    return(o)
  }

  PowerLink(
    link_name = paste0("power(lambda=", round(lambda, 5), ")"),
    link_bounds = c(0, Inf),
    link_params = list(lambda = lambda),
    lambda = lambda
  )
}
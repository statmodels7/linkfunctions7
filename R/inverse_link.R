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

# Exact analytical derivatives of the inverse link function (wrt eta)
# Due to the symmetric nature of f(x) = 1/x, these are structurally identical
# to the link function derivatives, but evaluated at eta.
S7::method(dlinkinv, InverseLink) <- function(x, eta) -1 / (eta^2)
S7::method(d2linkinv, InverseLink) <- function(x, eta)  2 / (eta^3)
S7::method(d3linkinv, InverseLink) <- function(x, eta) -6 / (eta^4)
S7::method(d4linkinv, InverseLink) <- function(x, eta) 24 / (eta^5)

#' @title The Inverse (Reciprocal) Link Function
#'
#' @include generics.R
#' @include link_class.R
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
#' @return An S7 object of class \code{InverseLink} (inheriting from \code{link}) containing the transformation functions
#' and their exact analytical derivatives up to the fourth order.
#'
#' @seealso \code{\link{link}}, \code{\link{identity_link}}
#' @export
inverse_link <- function() {
  InverseLink(
    link_name = "inverse",
    link_bounds = c(0, Inf),
    
    # The inverse link requires no additional mathematical parameters
    link_params = NULL
  )
}
#' @title S7 Class for the Probit Link
#'
#' @description
#' The class \code{\link{probit_link}} instantiates.
#'
#' @return An S7 object of class \code{ProbitLink}, inheriting from
#'   \code{\link{link}}.
#'
#' @seealso \code{\link{probit_link}}, the constructor users call.
#' @keywords internal
ProbitLink <- S7::new_class(
  name = "ProbitLink",
  parent = link
)

# --- Methods for ProbitLink ---

S7::method(linkfun, ProbitLink) <- function(x, theta) stats::qnorm(theta)
S7::method(linkinv, ProbitLink) <- function(x, eta) stats::pnorm(eta)

# Exact analytical derivatives of the link function (wrt theta)
# We calculate eta = qnorm(theta) and phi = dnorm(eta) locally 
# to significantly optimize computational operations.
S7::method(dlinkfun, ProbitLink) <- function(x, theta) lk_probit_fwd_cpp(theta, 1L)
S7::method(d2linkfun, ProbitLink) <- function(x, theta) lk_probit_fwd_cpp(theta, 2L)
S7::method(d3linkfun, ProbitLink) <- function(x, theta) lk_probit_fwd_cpp(theta, 3L)
S7::method(d4linkfun, ProbitLink) <- function(x, theta) lk_probit_fwd_cpp(theta, 4L)

# Exact analytical derivatives of the inverse link function (wrt eta)
# They rely exclusively on standard normal density properties.
S7::method(dlinkinv, ProbitLink) <- function(x, eta) lk_probit_inv_cpp(eta, 1L)
S7::method(d2linkinv, ProbitLink) <- function(x, eta) lk_probit_inv_cpp(eta, 2L)
S7::method(d3linkinv, ProbitLink) <- function(x, eta) lk_probit_inv_cpp(eta, 3L)
S7::method(d4linkinv, ProbitLink) <- function(x, eta) lk_probit_inv_cpp(eta, 4L)

#' @title The Probit Link Function
#'
#' @include generics.R
#' @include link_class.R
#' @description
#' The probit link \eqn{\eta = \Phi^{-1}(\theta)} on \eqn{(0, 1)}, with
#' \eqn{\Phi} the standard normal distribution function.
#' @details
#' The Probit link is mathematically defined as \eqn{\eta = \Phi^{-1}(\theta)}, where 
#' \eqn{\Phi^{-1}} is the quantile function of the standard normal distribution (\code{qnorm}).
#' The inverse link is \eqn{\theta = \Phi(\eta)}, the standard normal CDF (\code{pnorm}).
#'
#' Similarly to the \code{logit} link, the Probit is symmetric around \eqn{\theta = 0.5} 
#' (where \eqn{\eta = 0}). However, the tails of the Normal distribution approach 0 
#' and 1 faster than the Logistic distribution.
#'
#' The strictly mathematical domain of \eqn{\theta} is \code{c(0, 1)}.
#'
#' @return An S7 object of class \code{ProbitLink} (inheriting from \code{link}) containing the transformation functions
#' and their exact analytical derivatives up to the fourth order.
#'
#' @examples
#' lk <- probit_link()
#' lk
#'
#' p <- c(0.1, 0.5, 0.9)
#' eta <- linkfun(lk, p)      # standard normal quantiles
#' eta
#' linkinv(lk, eta)
#'
#' # the first inverse derivative is the standard normal density
#' dlinkinv(lk, 0)
#' dnorm(0)
#'
#' # probit tails approach 0 and 1 faster than logit ones
#' linkinv(probit_link(), 3)
#' linkinv(logit_link(), 3)
#'
#' @seealso \code{\link{link}}, \code{\link{logit_link}}, \code{\link{cauchit_link}}
#' @importFrom stats qnorm pnorm dnorm
#' @export
probit_link <- function() {
  ProbitLink(
    link_name = "probit",
    link_bounds = c(0, 1),
    link_params = NULL
  )
}

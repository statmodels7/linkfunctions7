#' @title S7 Class for the Power Link
#'
#' @description
#' Carries the power transformation \eqn{\eta = \theta^{\lambda}} on
#' \eqn{(0, \infty)} for a non-zero exponent, with inverse
#' \eqn{\theta = \eta^{1/\lambda}}. The exponent is stored in `link_params`, so
#' one class serves every \eqn{\lambda}.
#'
#' At \eqn{\lambda = 0} the constructor returns a [LogLink] instead, the limit
#' of the Box-Cox transformation \eqn{(\theta^\lambda - 1)/\lambda}.
#'
#' @param link_name A character string naming the link, set by the
#'   constructor and shown by `print()`.
#' @param link_bounds A length-two numeric vector, the open interval in
#'   which the parameter lies. Set by the constructor; see Value for this link's.
#' @param link_params A list of the link's own parameters, empty where it has
#'   none. Set by the constructor.
#' @param lambda The exponent of the transformation.
#'
#' @return An S7 object of class `PowerLink`, inheriting from [link()] and
#'   carrying its three properties `link_name`, `link_bounds` and
#'   `link_params`. Its `link_bounds` are `c(0, Inf)` and its `link_params` holds `lambda`.
#'
#' @section Methods:
#' Twelve methods are registered on this class: [linkfun()] and [linkinv()],
#' and the five derivative orders in each direction, [dlinkfun()] through
#' [d5linkfun()] going out and [dlinkinv()] through [d5linkinv()] coming
#' back. Both directions are a power, so all ten derivatives are falling
#' factorials in the exponent. Each is wrapped in [na_from()] because `NA^0`
#' is one in R, which would turn a missing parameter into a number as soon
#' as an exponent reached zero.
#'
#' @aliases linkfun.PowerLink
#' @aliases linkinv.PowerLink
#' @aliases dlinkfun.PowerLink
#' @aliases d2linkfun.PowerLink
#' @aliases d3linkfun.PowerLink
#' @aliases d4linkfun.PowerLink
#' @aliases d5linkfun.PowerLink
#' @aliases dlinkinv.PowerLink
#' @aliases d2linkinv.PowerLink
#' @aliases d3linkinv.PowerLink
#' @aliases d4linkinv.PowerLink
#' @aliases d5linkinv.PowerLink
#'
#' @seealso [power_link()], the constructor users call.
#' @keywords internal
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

S7::method(dlinkfun, PowerLink) <- function(x, theta) na_from(x@lambda * (theta^(x@lambda - 1)), theta)
S7::method(d2linkfun, PowerLink) <- function(x, theta) na_from(x@lambda * (x@lambda - 1) * (theta^(x@lambda - 2)), theta)
S7::method(d3linkfun, PowerLink) <- function(x, theta) na_from(x@lambda * (x@lambda - 1) * (x@lambda - 2) * (theta^(x@lambda - 3)), theta)
S7::method(d4linkfun, PowerLink) <- function(x, theta) na_from(x@lambda * (x@lambda - 1) * (x@lambda - 2) * (x@lambda - 3) * (theta^(x@lambda - 4)), theta)
S7::method(d5linkfun, PowerLink) <- function(x, theta) na_from(x@lambda * (x@lambda - 1) * (x@lambda - 2) * (x@lambda - 3) * (x@lambda - 4) * (theta^(x@lambda - 5)), theta)

S7::method(dlinkinv, PowerLink) <- function(x, eta) { k <- 1 / x@lambda; na_from(k * (eta^(k - 1)), eta) }
S7::method(d2linkinv, PowerLink) <- function(x, eta) { k <- 1 / x@lambda; na_from(k * (k - 1) * (eta^(k - 2)), eta) }
S7::method(d3linkinv, PowerLink) <- function(x, eta) { k <- 1 / x@lambda; na_from(k * (k - 1) * (k - 2) * (eta^(k - 3)), eta) }
S7::method(d4linkinv, PowerLink) <- function(x, eta) { k <- 1 / x@lambda; na_from(k * (k - 1) * (k - 2) * (k - 3) * (eta^(k - 4)), eta) }
S7::method(d5linkinv, PowerLink) <- function(x, eta) { k <- 1 / x@lambda; na_from(k * (k - 1) * (k - 2) * (k - 3) * (k - 4) * (eta^(k - 5)), eta) }

#' @title The Power Link Function
#'
#' @include generics.R
#' @include link_class.R
#' @description
#' The power link \eqn{\eta = \theta^\lambda} on \eqn{(0, \infty)}; at
#' \eqn{\lambda = 0} it returns the log link, the limit of the Box-Cox
#' transformation.
#' @param lambda A numeric value defining the power of the transformation. Defaults to 1.
#'
#' @details
#' The power link is \eqn{\eta = \theta^\lambda}, with inverse
#' \eqn{\theta = \eta^{1/\lambda}}.
#'
#' As \eqn{\lambda \to 0} the Box-Cox transformation
#' \eqn{(\theta^\lambda - 1)/\lambda} tends to \eqn{\log(\theta)}, so for
#' `lambda = 0` the function returns a [log_link()] object whose
#' `link_params` record `lambda = 0`.
#'
#' Common special cases are
#'
#' - `lambda = 1`: identity link;
#' - `lambda = 0.5`: square-root link;
#' - `lambda = -1`: inverse link;
#' - `lambda = 0`: log link.
#'
#' The domain of \eqn{\theta} is \eqn{(0, \infty)}. For every non-zero
#' `lambda` the image of the link is also \eqn{(0, \infty)}, so the linear
#' predictor must stay positive during optimization; [eta_bounds()] returns
#' this range.
#'
#' @return An S7 object of class `PowerLink`, inheriting from [link()], or of
#'   class `LogLink` when `lambda = 0`.
#'
#' @examples
#' lk <- power_link(2)
#' lk
#'
#' theta <- c(1, 2, 3)
#' eta <- linkfun(lk, theta)
#' eta
#' linkinv(lk, eta)
#'
#' # special cases of the family
#' linkfun(power_link(1),    5)   # identity
#' linkfun(power_link(0.5),  4)   # square root
#' linkfun(power_link(-1),   4)   # inverse
#'
#' # lambda = 0 returns the log link, the limit of the Box-Cox transformation
#' power_link(0)
#' linkfun(power_link(0), exp(1))
#'
#' @seealso [link()], [log_link()], [identity_link()]
#' @export
power_link <- function(lambda = 1) {
  # Without this the `if` below is what fails, on a missing or vector lambda,
  # with a message about the condition rather than about the argument.
  if (length(lambda) != 1L || !is.numeric(lambda) || !is.finite(lambda)) {
    stop("'lambda' must be a single finite number.", call. = FALSE)
  }

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

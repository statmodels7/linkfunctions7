#' @include link_class.R bounded_link.R power_link.R softplus_link.R
NULL

#' The Scalar Route of a Link
#'
#' @description
#' Returns how the package's compiled scalar entry points address a link: the
#' class name that the C function `lf7_class_id()` recognizes, and the link's
#' own parameters in the order `lf7_inv12p()` reads them. A consumer that
#' resolves the entry points with `R_GetCCallable()` uses the name to obtain
#' the link's identifier once and passes the parameters at every evaluation.
#'
#' @details
#' The entry points are `lf7_class_id(name)` and `lf7_inv12p(id, par, eta,
#' h, h1, h2)`, which writes the inverse link and its first two derivatives
#' at one value of the linear predictor, each as the class's R method
#' computes it, and `lf7_clamp(theta, lwr, upr)`, which applies
#' [link_bounds_clamp()] to one value. The parameters are \eqn{\lambda} for
#' [power_link()], \eqn{a} for [softplus_link()], the lower bound and the
#' width for a doubly bounded [bounded_link()], and the single bound for a
#' one-sided one; the other links have none.
#'
#' The class is matched by its exact name, so a class defined outside the
#' package, including one that inherits from a class of the package and may
#' override its methods, returns `NULL`.
#'
#' @param x A link object inheriting from [link()].
#'
#' @return `NULL` when the entry points do not cover the link; otherwise a
#'   list with `name`, a single character string, and `par`, a numeric
#'   vector of the link's parameters (of length zero for a link without
#'   parameters).
#'
#' @examples
#' link_scalar_route(logit_link())
#' link_scalar_route(bounded_link(lwr = 2, upr = 5))
#' link_scalar_route(power_link(0.5))
#'
#' @export
link_scalar_route <- function(x) {
  name <- attr(S7::S7_class(x), "name")
  par <- switch(name,
    IdentityLink = , LogLink = , LogitLink = , ProbitLink = ,
    ClogLogLink = , LogLogLink = , CauchitLink = , RhobitLink = ,
    SqrtLink = , InverseLink = , InverseSqLink = numeric(0),
    PowerLink = x@lambda,
    SoftplusLink = x@a,
    DoublyBoundedLink = c(x@lwr, x@width),
    LowerBoundedLink = x@lwr,
    UpperBoundedLink = x@upr,
    NULL)
  if (is.null(par)) return(NULL)
  list(name = name, par = as.numeric(par))
}

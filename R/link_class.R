#' @title S7 Class for Statistical Link Functions
#'
#' @import S7
#' @description
#' The base S7 class for link functions. It carries the name, the domain and
#' any link parameters; the transformations themselves -- forward, inverse and
#' their analytical derivatives to fourth order -- are methods that each
#' subclass registers on the ten generics.
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
#'
#' @return An S7 object of class \code{link}. In practice this class is not
#'   instantiated directly: each link is a subclass created by one of the
#'   constructors (\code{\link{logit_link}}, \code{\link{power_link}}, ...), and
#'   \code{link} is what they all inherit from and what methods dispatch on.
#'
#' @examples
#' # every constructor returns an object inheriting from `link`
#' lk <- logit_link()
#' lk
#' S7::S7_inherits(lk, link)
#'
#' lk@link_name
#' lk@link_bounds
#'
#' @export
link <- S7::new_class(
  name = "link",
  properties = list(
    link_name = S7::class_character,
    link_bounds = S7::class_numeric,
    link_params = S7::class_any
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


#' A Constant Vector That Preserves Missingness
#'
#' @description
#' Returns \code{value} repeated to the length of \code{v}, but missing wherever
#' \code{v} is missing.
#'
#' @details
#' A derivative that reduces to a constant must still report that it does not
#' know the answer for an input it was not given. R makes this easy to get wrong:
#' \code{NA^0} is \code{1}, so \code{theta^(lambda - 2)} silently turns a missing
#' parameter into a number as soon as \code{lambda} is 2. Every derivative method
#' that returns a constant (the identity link's, the square root link's third and
#' fourth inverse derivatives) goes through this helper instead of \code{rep()}.
#'
#' @param v A numeric vector whose length and missingness pattern are copied.
#' @param value The constant to repeat.
#'
#' @return A numeric vector as long as \code{v}, equal to \code{value} except
#'   where \code{v} is \code{NA}.
#'
#' @seealso \code{\link{na_from}}, the same idea for a computed result.
#' @keywords internal
const_like <- function(v, value) {
  out <- rep(value, length(v))
  out[is.na(v)] <- NA_real_
  out
}

#' Carry Missingness From an Input Over to a Result
#'
#' @description
#' Sets \code{r} to \code{NA} wherever \code{v} is \code{NA}.
#'
#' @details
#' Same hazard as \code{\link{const_like}}, one step further along: an expression
#' whose exponent happens to vanish stops depending on its argument, and loses the
#' argument's missingness along with it. The power link is the affected case,
#' \code{theta^(lambda - 2)} being exactly \code{1} for a missing \code{theta}
#' once \code{lambda} is 2.
#'
#' @param r A numeric vector, the computed result.
#' @param v The numeric vector the result was computed from.
#'
#' @return \code{r}, with \code{NA} in every position where \code{v} is \code{NA}.
#'
#' @seealso \code{\link{const_like}}
#' @keywords internal
na_from <- function(r, v) {
  r[is.na(v)] <- NA_real_
  r
}

#' Derivatives of the Standard Logistic Function
#'
#' @description
#' The \code{k}-th derivative of \eqn{\sigma(z) = 1/(1 + e^{-z})}, written as a
#' polynomial in \eqn{p = \sigma(z)} itself.
#'
#' @details
#' Three separate links need these same four polynomials, which is why they live
#' in one place rather than being written out three times:
#' \itemize{
#'   \item \code{\link{logit_link}} uses them directly, \eqn{h^{(k)} = \sigma^{(k)}};
#'   \item \code{\link{bounded_link}} with both endpoints scales them by the
#'     interval width, \eqn{h^{(k)} = W \sigma^{(k)}};
#'   \item \code{\link{softplus_link}} uses them shifted one order down, since the
#'     softplus is an antiderivative of the logistic: \eqn{h^{(k+1)} = a^k \sigma^{(k)}}.
#' }
#' The polynomials are
#' \deqn{\sigma' = p(1-p)}
#' \deqn{\sigma'' = p(1-p)(1-2p)}
#' \deqn{\sigma''' = p(1-p)(1 - 6p + 6p^2)}
#' \deqn{\sigma'''' = p(1-p)(1 - 14p + 36p^2 - 24p^3)}
#' and are evaluated in Horner form, which is twice as fast at the fourth order
#' and agrees with the expanded form to within one unit in the last place.
#'
#' @param p A numeric vector of logistic values, \eqn{p = \sigma(z)}.
#' @param k The derivative order, an integer from 1 to 4.
#'
#' @return A numeric vector of the same length as \code{p}.
#'
#' @keywords internal
logistic_deriv <- function(p, k) {
  pq <- p * (1 - p)
  switch(k,
    pq,
    pq * (1 - 2 * p),
    pq * (1 + p * (-6 + 6 * p)),
    pq * (1 + p * (-14 + p * (36 - 24 * p)))
  )
}

#' The Smallest Parameter Value the Exponential Links Will Report
#'
#' @description
#' The floor applied to \code{exp(eta)} by every link whose inverse is an
#' exponential (\code{\link{log_link}}, \code{\link{cloglog_link}}, and the
#' lower- and upper-bounded links).
#'
#' @details
#' The floor exists so that a parameter reported as \eqn{\theta} can be divided
#' into without producing \code{Inf}: the forward derivatives of these links are
#' \eqn{1/\theta}, \eqn{-1/\theta^2}, \eqn{2/\theta^3} and \eqn{-6/\theta^4}, and
#' the fourth is the binding one. Solving \eqn{6/\theta^4 \le} \code{double.xmax}
#' and keeping a factor of four in hand gives
#' \code{(24 / .Machine$double.xmax)^0.25}, about \code{1.9e-77}, at which
#' \eqn{-6/\theta^4} evaluates to \code{-4.5e307}.
#'
#' The point of choosing it this way is that the floor should be as \emph{low} as
#' that constraint allows, not as high as seems safe. It was previously
#' \code{.Machine$double.eps}, which is 61 orders of magnitude higher than
#' necessary and silently corrupted \eqn{\theta} for every \eqn{\eta < -36}:
#' \code{linkinv(log_link(), -40)} returned \code{2.2e-16} instead of
#' \code{4.2e-18}, and the round trip came back \code{-36.04} instead of
#' \code{-40}. The present value keeps \eqn{\theta} exact down to
#' \eqn{\eta \approx -177} while leaving every derivative just as finite as
#' before.
#'
#' @format A length-one numeric vector.
#' @return A length-one numeric vector, about \code{1.9e-77}. This is a constant
#'   and not a function, but a documented topic needs a return section either
#'   way: it is one of the two things a first CRAN submission is most often sent
#'   back for.
#' @seealso \code{\link{exp_floored}}
#' @keywords internal
exp_floor <- (24 / .Machine$double.xmax)^0.25

#' A Floored Exponential
#'
#' @description
#' \code{exp(eta)}, bounded below by \code{\link{exp_floor}}.
#'
#' @param eta A numeric vector of linear predictors.
#'
#' @return A numeric vector, never smaller than \code{\link{exp_floor}}.
#'
#' @keywords internal
exp_floored <- function(eta) pmax(exp(eta), exp_floor)


#' Clamp a Parameter Strictly Inside Its Domain
#'
#' @description
#' Moves a value that has reached or passed a bound to the nearest double
#' strictly inside it, and a non-finite value to the largest finite double of
#' that sign. Applied by \code{\link{linkinv}} to every link.
#'
#' @details
#' A link is documented as a bijection onto an \strong{open} interval, and in
#' exact arithmetic it is. In double precision it is not: \code{plogis(37)} is
#' exactly 1, \code{2 + exp(-40)} is exactly 2, and \code{exp(800)} is infinite.
#' A caller then receives a probability of exactly 1, or a variance of exactly
#' 0, and the next thing they do is take its logarithm or divide by it.
#'
#' The correction is the smallest one that can work and it is derived, not
#' chosen. What binds is "strictly inside, and finite",
#' and its two extremes are the neighboring representable double and the
#' largest finite one.
#'
#' \subsection{The relative bump}{
#' R has no \code{nextafter}, and the arithmetic substitute has to respect that
#' \strong{the spacing of doubles is absolute near a non-zero bound}. One ulp at
#' 2 is about 4.4e-16 while one ulp at 1e-300 is about 1e-316, so a single
#' additive constant cannot serve both. \code{b + |b| * eps} is one to two ulps
#' from \code{b} at any magnitude, which is strictly inside and as close as
#' arithmetic reliably gets.
#'
#' A bound at zero is the exception and needs no bump, since there the spacing
#' is relative all the way down to 1e-308 and the exponential links already
#' floor at \code{\link{exp_floor}}. The clamp therefore leaves an exact zero
#' bound to \code{\link{exp_floor}} and uses the smallest positive normal only
#' if something has still landed on it.
#' }
#'
#' @param theta A numeric vector, as a method computed it.
#' @param bounds The link's \code{link_bounds}, a length-2 numeric vector.
#'
#' @return \code{theta}, with any value that has landed exactly on a bound
#'   moved just inside it and any infinity brought back to the largest finite
#'   double. \code{NA} and \code{NaN} pass through untouched, as does a value
#'   strictly outside by a real margin: converting either would hide something
#'   rather than fix it.
#'
#' @examples
#' link_bounds_clamp(c(0, 0.5, 1), c(0, 1))
#' link_bounds_clamp(c(2, 3, Inf), c(2, Inf))
#'
#' @seealso \code{\link{linkinv}}
#' @export
link_bounds_clamp <- function(theta, bounds) {
  lwr <- bounds[1]
  upr <- bounds[2]
  eps <- .Machine$double.eps
  big <- .Machine$double.xmax

  # Infinities first, so that the comparisons below see a number. NaN is left
  # alone deliberately, and Inf is a value a caller cannot use either way.
  inf <- is.infinite(theta)
  if (any(inf)) theta[inf] <- sign(theta[inf]) * big

  # EXACTLY on the bound, not merely outside it. Saturation is the arithmetic
  # running out of resolution while eta was a perfectly good input, and it lands
  # the result precisely on the endpoint every time: plogis(37) is 1, not 1 plus
  # something. A value strictly outside by a real margin is a different
  # question -- inverse_link() at eta = -40 returns -0.025, because 1/eta is a
  # bijection from (0, Inf) and -40 is not an admissible linear predictor for
  # it. Clamping that would turn "you gave me an eta this link cannot take" into
  # a small positive number, which is worse than the complaint.
  if (is.finite(lwr)) {
    low <- !is.na(theta) & theta == lwr
    if (any(low)) {
      theta[low] <- if (lwr == 0) .Machine$double.xmin else lwr + abs(lwr) * eps
    }
  }
  if (is.finite(upr)) {
    high <- !is.na(theta) & theta == upr
    if (any(high)) {
      theta[high] <- if (upr == 0) -.Machine$double.xmin else upr - abs(upr) * eps
    }
  }
  theta
}

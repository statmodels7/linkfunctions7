#' @title S7 Class for Statistical Link Functions
#'
#' @import S7
#' @description
#' The base S7 class for link functions. It carries the name, the domain and
#' any link parameters. The transformations themselves are methods that each
#' subclass registers on the generics for the forward map, the inverse, and
#' their analytical derivatives to fifth order in both directions.
#'
#' @details
#' Objects of class `link` are instantiated using the S7 object system.
#'
#' The documentation uses the following notation:
#'
#' - \eqn{\theta}: The response parameter (e.g., probability, mean, dispersion).
#' - \eqn{\eta}: The linear predictor (unconstrained scale).
#'
#' The relationship is defined as \eqn{\eta = g(\theta)} (link function) and
#' \eqn{\theta = g^{-1}(\eta)} (inverse link function).
#'
#' @param link_name A character string identifying the link (e.g., "logit").
#' @param link_bounds A numeric vector of length 2 `c(lower, upper)` defining the valid domain for \eqn{\theta}.
#' @param link_params A list or vector of additional parameters required to define the link, or `NULL`.
#'
#' @return An S7 object of class `link`. In practice this class is not
#'   instantiated directly: each link is a subclass created by one of the
#'   constructors ([logit_link()], [power_link()], ...), and
#'   `link` is what they all inherit from and what methods dispatch on.
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
#' @seealso [linkfun()], [linkinv()], [linkderiv()], [linkinvderiv()], [check_link()]
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
#' Returns `value` repeated to the length of `v`, but missing wherever
#' `v` is missing.
#'
#' @details
#' A derivative that reduces to a constant must still return `NA` where its
#' input is `NA`. This is easy to get wrong in R: `NA^0` is `1`, so
#' `theta^(lambda - 2)` silently turns a missing parameter into a number as
#' soon as `lambda` is 2. Every derivative method that returns a constant (the
#' identity link's, and the square root link's second to fifth inverse
#' derivatives) goes through this helper instead of `rep()`.
#'
#' @param v A numeric vector whose length and missingness pattern are copied.
#' @param value The constant to repeat.
#'
#' @return A numeric vector as long as `v`, equal to `value` except
#'   where `v` is `NA`.
#'
#' @seealso [na_from()], the same idea for a computed result.
#' @keywords internal
const_like <- function(v, value) {
  out <- rep(value, length(v))
  out[is.na(v)] <- NA_real_
  out
}

#' Carry Missingness From an Input Over to a Result
#'
#' @description
#' Sets `r` to `NA` wherever `v` is `NA`.
#'
#' @details
#' This handles the same problem as [const_like()] for a computed result: an
#' expression whose exponent vanishes no longer depends on its argument, and
#' loses the argument's missingness with it. The power link is the case in
#' point, since `theta^(lambda - 2)` is exactly `1` for a missing `theta`
#' once `lambda` is 2.
#'
#' @param r A numeric vector, the computed result.
#' @param v The numeric vector the result was computed from.
#'
#' @return `r`, with `NA` in every position where `v` is `NA`.
#'
#' @seealso [const_like()]
#' @keywords internal
na_from <- function(r, v) {
  r[is.na(v)] <- NA_real_
  r
}

#' Derivatives of the Standard Logistic Function
#'
#' @description
#' The `k`-th derivative of \eqn{\sigma(z) = 1/(1 + e^{-z})}, written as a
#' polynomial in \eqn{p = \sigma(z)} itself.
#'
#' @details
#' Three links use these polynomials:
#'
#' - [logit_link()] uses them directly, \eqn{h^{(k)} = \sigma^{(k)}};
#' - [bounded_link()] with both endpoints scales them by the
#'   interval width, \eqn{h^{(k)} = W \sigma^{(k)}};
#' - [softplus_link()] uses them shifted one order down, since the
#'   softplus is an antiderivative of the logistic: \eqn{h^{(k+1)} = a^k \sigma^{(k)}}.
#'
#' The three links call the transcription of these polynomials in
#' `src/link_kernels.cpp`, which replaced the R bodies when the
#' transcendental links were compiled. This function is the R statement of
#' the same five polynomials. `test-logistic-twin.R` compares the two at
#' every order and checks that each of the three links reaches the
#' polynomial named in its description, and `lk_logistic_poly_cpp()` calls
#' the compiled version directly.
#'
#' The comparison uses a tolerance instead of exact equality, because both
#' forms are Horner evaluations that contain multiply-adds, and a compiler may
#' fuse these into FMA instructions that skip an intermediate rounding.
#'
#' The polynomials are
#' \deqn{\sigma' = p(1-p)}
#' \deqn{\sigma'' = p(1-p)(1-2p)}
#' \deqn{\sigma''' = p(1-p)(1 - 6p + 6p^2)}
#' \deqn{\sigma'''' = p(1-p)(1 - 14p + 36p^2 - 24p^3)}
#' \deqn{\sigma^{(5)} = p(1-p)(1 - 30p + 150p^2 - 240p^3 + 120p^4)}
#' and are evaluated in Horner form, which is faster than the expanded form
#' and agrees with it up to rounding error; near a root of a polynomial that
#' error can be large relative to the value. Each polynomial
#' follows from the one before by \eqn{P_{k+1} = (1-2p)P_k + p(1-p)P_k'}, so a
#' higher order can be generated from the recurrence.
#'
#' @param p A numeric vector of logistic values, \eqn{p = \sigma(z)}.
#' @param k The derivative order, an integer from 1 to 5.
#'
#' @return A numeric vector of the same length as `p`.
#'
#' @keywords internal
logistic_deriv <- function(p, k) {
  pq <- p * (1 - p)
  switch(k,
    pq,
    pq * (1 - 2 * p),
    pq * (1 + p * (-6 + 6 * p)),
    pq * (1 + p * (-14 + p * (36 - 24 * p))),
    pq * (1 + p * (-30 + p * (150 + p * (-240 + 120 * p))))
  )
}

#' The Floor of the Exponential Links
#'
#' @description
#' The floor applied to `exp(eta)` by every link whose inverse is an
#' exponential ([log_link()], [cloglog_link()], and the
#' lower- and upper-bounded links).
#'
#' @details
#' The floor exists so that a parameter returned as \eqn{\theta} can be
#' divided into without producing `Inf`: the forward derivatives of these
#' links are \eqn{1/\theta}, \eqn{-1/\theta^2}, \eqn{2/\theta^3} and
#' \eqn{-6/\theta^4}, and the fourth is the binding one. Solving
#' \eqn{6/\theta^4 \le} `double.xmax` and keeping a factor of four in hand
#' gives `(24 / .Machine$double.xmax)^0.25`, about `1.9e-77`, at which
#' \eqn{-6/\theta^4} evaluates to `-4.5e307`.
#'
#' The floor is the lowest value that this constraint allows, and it keeps
#' \eqn{\theta} exact down to \eqn{\eta \approx -176.7}. The fifth forward
#' derivative, \eqn{24/\theta^5}, does not fit under it and is `Inf` for
#' \eqn{\eta} below about \eqn{-141}.
#'
#' @format A length-one numeric vector.
#' @return A length-one numeric vector, about `1.9e-77`, at which
#'   \eqn{-6/\theta^4} evaluates to `-4.5e307`.
#' @seealso [exp_floored()]
#' @keywords internal
exp_floor <- (24 / .Machine$double.xmax)^0.25

#' A Floored Exponential
#'
#' @description
#' `exp(eta)`, bounded below by [exp_floor()].
#'
#' @param eta A numeric vector of linear predictors.
#'
#' @return A numeric vector, never smaller than [exp_floor()].
#'
#' @keywords internal
exp_floored <- function(eta) pmax(exp(eta), exp_floor)


#' Clamp a Parameter Strictly Inside Its Domain
#'
#' @description
#' Moves a value that lies exactly on a finite bound to a double strictly
#' inside it, one or two units in the last place away, and an infinite value
#' to the largest finite double of that sign. Applied by [linkinv()] to every link.
#'
#' @details
#' In exact arithmetic a link is a bijection onto an open interval, but in
#' double precision the endpoints can be reached: `plogis(37)` is exactly 1,
#' `2 + exp(-40)` is exactly 2, and `exp(800)` is infinite. A caller would then
#' receive a probability of exactly 1 or a variance of exactly 0, and taking
#' its logarithm or dividing by it fails.
#'
#' The clamp moves such a value by the smallest amount that makes it strictly
#' inside the interval and finite: to a neighboring representable double, or
#' to the largest finite double.
#'
#' \subsection{The size of the step}{
#' R has no `nextafter()`, so the step is computed arithmetically. Near a
#' non-zero bound \eqn{b} the spacing of doubles is set by the magnitude of
#' \eqn{b}: one ulp at 2 is about 4.4e-16, while one ulp at 1e-300 is about
#' 1e-316, so no single additive constant serves both. The value
#' `b + |b| * eps` lies one to two ulps from `b` at any magnitude, which is
#' strictly inside the interval and as close to `b` as the arithmetic reliably
#' allows.
#'
#' A bound at zero needs no such step, because the exponential links already
#' floor their result at [exp_floor()]. A value that still lands exactly on a
#' zero bound is moved to the smallest positive normal double.
#' }
#'
#' @param theta A numeric vector, as a method computed it.
#' @param bounds The link's `link_bounds`, a length-2 numeric vector.
#'
#' @return `theta`, with any value that has landed exactly on a bound
#'   moved just inside it and any infinity brought back to the largest finite
#'   double. `NA` and `NaN` pass through unchanged, as does a value outside
#'   the bounds by more than rounding, since changing either would hide an
#'   error in the calling code.
#'
#' @examples
#' link_bounds_clamp(c(0, 0.5, 1), c(0, 1))
#' link_bounds_clamp(c(2, 3, Inf), c(2, Inf))
#'
#' @seealso [linkinv()]
#' @export
link_bounds_clamp <- function(theta, bounds) {
  lwr <- bounds[1]
  upr <- bounds[2]
  eps <- .Machine$double.eps
  big <- .Machine$double.xmax

  # This body runs from linkinv()'s generic for EVERY link on every call, so
  # it looks like the place to spend an optimization, and a range() over
  # theta does decide all four questions below at once without allocating.
  # Measured, it is not worth doing: on a two-sided bound it gains about
  # 1.6x on the call, on a one-sided bound it loses, and end to end a
  # gaussian fit at n = 100000 got 40 per cent SLOWER, the extra pass
  # costing more than the logical vectors it avoided. The elementwise form
  # is kept.

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

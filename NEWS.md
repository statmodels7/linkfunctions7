# linkfunctions7 0.1.0

## Numerical behavior

* The derivatives of the transcendental links -- logit, probit, cloglog,
  loglog, cauchit, rhobit and softplus -- are computed by compiled kernels
  in `src/link_kernels.cpp`. The doubly bounded link and the softplus
  inverse reuse the logit kernel, scaled by the width and shifted one
  order. Measured through S7 dispatch, `d4linkinv()` on the logit falls
  from 40 ms to 25 ms at a million points; what the port buys is not the
  transcendental, which costs the same in either language, but the vector
  temporaries an order-four polynomial allocates.

* The finite-difference weights and offsets come from `numericals7`, where
  the toolkit keeps one stencil library.

* `linkinv()` returns a value strictly inside the open parameter interval.
  Nine of the eighteen links previously reached a bound or went non-finite
  under a large predictor: `linkinv(logit_link(), 37)` was exactly 1 and
  the round trip through `linkfun()` then `Inf`, which `distributions7`
  rejects, since it validates against open intervals. The correction is
  made in the generic, so every link inherits it and a user-written link
  needs no change. It fires only on exact equality with a bound; a value
  outside by a real margin is left alone, that being a complaint rather
  than a rounding artifact.

* `cloglog_link()`'s forward direction is computed through `log1p`. Written
  as `log(-log(1 - theta))` it returned `-Inf` for small `theta` where the
  value is finite and representable (-176.66 at `theta = 1.9e-77`), taking
  the four forward derivatives that divide by that logarithm with it.

* The floor the exponential links apply is derived from the quantity that
  binds rather than taken from a familiar constant. The old value,
  `.Machine$double.eps`, is the resolution of a double near one and not a
  lower bound on anything relevant; it corrupted `theta` for every
  `eta < -36`. The constraint is that the log link's fourth derivative
  `-6/theta^4` must not overflow, giving `(24/double.xmax)^(1/4)`, about
  1.9e-77. The same correction applies to `cloglog` and the lower-bounded
  link.

* `softplus_link()` is written in `u = -expm1(-a*theta)` rather than in
  `expm1(a*theta)`. The old form gave `NaN` for `d4linkfun` from
  `theta = 177` at `a = 1` and from `theta = 24` at `a = 30`, an ordinary
  value for a steep softplus.

## Derivatives

* A link that implements only `linkfun()` and `linkinv()` is complete: the
  base class supplies the eight derivative generics numerically. Each
  fallback finds the highest order implemented analytically and applies one
  central stencil of order `k - m` to it, never a chain of first
  differences. On a log link defined with nothing the fourth derivative is
  accurate to about 2e-5; supplying the analytic first and second takes it
  to 2e-8.

* S7 classes are compared by name and package rather than by object
  identity. `identical()` on a class is `FALSE` for a class re-created from
  the same definition, which is what happens whenever the code is
  re-evaluated rather than loaded, so the base numerical fallback was taken
  for an analytic method and every fallback differentiated the order below
  it: the log link's fourth derivative was wrong by a factor of 900.

## Validation and documentation

* `check_link()` reports an unimplemented derivative as a failure. An
  order that raises left `NA` and broke the loop, and the summary reduced
  with `na.rm = TRUE`, so a link implementing only its first derivative was
  reported as passing all four orders.

* `check_link()` distinguishes an order it verified from one it could not.
  An order supplied by a numerical fallback would be compared against a
  numerical differentiation of the order below, which is the same
  arithmetic twice and agrees however wrong the link is; such an order is
  left `NA` and printed as `[numerical]`. `link_fallback_orders()` answers
  the same question programmatically.

* Every object in the namespace has a help page, every exported topic a
  `\value` and an executable example, and `tests/testthat/test-docs.R`
  asks all of that plus whether the pkgdown index is complete.

* A vignette on working with link functions, a README with badges, a
  pkgdown site and continuous integration on five platforms.

# Changelog

## linkfunctions7 0.5.0

- [`link_scalar_route()`](https://statmodels7.github.io/linkfunctions7/reference/link_scalar_route.md)
  returns the name and the parameters by which the compiled registry
  knows a link, or `NULL`. Every link class has scalar C entries,
  `lf7_class_id()` and `lf7_inv12p()` beside `lf7_scalar_id()`,
  registered with `R_RegisterCCallable()`, which a compiled consumer
  such as the `gas()` filter of modelterms7 reads instead of calling
  back into R. The vector kernels and the entries share
  `src/link_points.h`.

## linkfunctions7 0.4.0

- Every link carries its **fifth** derivative analytically, in both
  directions:
  [`d5linkfun()`](https://statmodels7.github.io/linkfunctions7/reference/d5linkfun.md)
  and
  [`d5linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/d5linkinv.md),
  registered on all sixteen classes, with
  [`linkderiv()`](https://statmodels7.github.io/linkfunctions7/reference/linkderiv.md)
  and
  [`linkinvderiv()`](https://statmodels7.github.io/linkfunctions7/reference/linkinvderiv.md)
  routing to order 5 and
  [`check_link()`](https://statmodels7.github.io/linkfunctions7/reference/check_link.md)
  validating it. The order exists because each order of differentiating
  a score-driven filter’s predictor through its own recursion draws in
  one more order of the link and of the family – the curvature reaches
  the third, the directional third derivative the fourth, and the outer
  Hessian of a model carrying such a term the fifth.

  Nothing was transcribed on trust. Each formula was derived from the
  order below it and checked against one Richardson pass on the analytic
  fourth before it was written, with a negative control on the same
  grid: a coefficient 5 per cent out reads between 7.9e-03 and 3.9 where
  the formula as shipped reads between 4.5e-12 and 7.4e-11. The
  recurrences are worth keeping, since they generate the orders above
  rather than being copied: the logistic polynomials follow , the
  probit’s inverse derivatives are , the cloglog’s obey and the
  softplus’s forward numerators are the Eulerian triangle, 1, 11, 11, 1
  at this order.

  The softplus and the doubly bounded link have no algebra of their own
  here: the first uses the logistic polynomials one order down times a
  power of its steepness, the second scales them by the interval width,
  and both identities are asserted rather than assumed.

- An unsupported order asks a compiled kernel for a number and gets
  `NA`. Every kernel takes the order as an argument and reached its
  order-4 branch through `default:`, so `lk_logit_inv_cpp(eta, 5L)` was
  [`identical()`](https://rdrr.io/r/base/identical.html) to the fourth
  order and `lk_probit_inv_cpp(eta, 9L)` answered as well. The routers
  have always rejected an order they do not carry, so nothing was
  reachable from the public surface and no result was ever wrong; the
  shape was, and it would have become live with the fifth order wired
  in.

- ⚠️ **The fifth forward derivative of the exponential links does not
  fit under the exponential floor, and the floor stays where it was
  derived.** `exp_floor` is , chosen so that cannot overflow, and at the
  floor that derivative reads -4.49e+307 against a `double.xmax` of
  1.80e+308: the bound is tight. The fifth would need = 5.82e-62,
  fifteen orders higher, so
  [`d5linkfun()`](https://statmodels7.github.io/linkfunctions7/reference/d5linkfun.md)
  of a log, lower- or upper-bounded link is `Inf` below about . It is
  not a cancellation and there is no rewrite: at the floor the value of
  is about 3.6e+384, which is not a double. Raising the floor to suit
  would clamp for between -177 and -141, where
  [`linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/linkinv.md)
  is exact today. The inverse direction, which is the one the chain rule
  onto the link scale uses, stays finite at every order, and a test pins
  both halves.

## linkfunctions7 0.3.0

- [`logistic_deriv()`](https://statmodels7.github.io/linkfunctions7/reference/logistic_deriv.md)
  is held to the compiled kernel it is the R statement of. The four
  logistic derivative polynomials are written twice, in R and as the
  `static inline logistic_poly()` in `src/link_kernels.cpp` that the
  logit, the doubly bounded link and the softplus reach; nothing called
  the R one and no test compared them, so the page’s claim that the two
  agree was unheld. `test-logistic-twin.R` compares them at all four
  orders over 501 points and checks that each of the three links reaches
  the polynomial its description names – as they stand, scaled by the
  interval width, and one order down times a power of the steepness.

  The comparison carries a tolerance rather than asking for identity.
  Both forms are Horner and contain multiply-adds, which a compiler may
  contract into an FMA; measured here the two agree to the bit at every
  order, and that is a property of one compiler. A negative control
  asserts that a polynomial wrong in one coefficient fails by more than
  1e-3, which no contraction reaches.

- [`eta_bounds()`](https://statmodels7.github.io/linkfunctions7/reference/eta_bounds.md)
  is exported. It answers what a link maps **from**, where `link_bounds`
  says what it maps **onto**, and the two are different questions:
  [`sqrt_link()`](https://statmodels7.github.io/linkfunctions7/reference/sqrt_link.md),
  [`inverse_link()`](https://statmodels7.github.io/linkfunctions7/reference/inverse_link.md),
  [`inverse_sq_link()`](https://statmodels7.github.io/linkfunctions7/reference/inverse_sq_link.md)
  and
  [`power_link()`](https://statmodels7.github.io/linkfunctions7/reference/power_link.md)
  at a positive exponent all reach the positive half of the theta axis
  from the positive half of the eta axis alone. A consumer that carries
  an unconstrained vector and applies an inverse link to each coordinate
  needs the second question answered, the map being neither defined nor
  injective outside those bounds – `linkinv(sqrt_link(), -2)` and
  `linkinv(sqrt_link(), 2)` are both 4 and the round trip returns the
  absolute value.

  Nothing about the function changes. It was already documented and
  reachable through `:::`; what changes is that a package depending on
  this one can now ask the question at construction, where the message
  can name the link, instead of discovering the fold in a fitted number.

## linkfunctions7 0.2.0

- Scalar C entry points for the fast route of a score-driven filter
  (piano_parallel.txt, section 2a), registered with
  `R_RegisterCCallable`: `lf7_scalar_id` (identity and log; an unknown
  name answers -1 and the consumer keeps its R callbacks), `lf7_inv12`
  (the inverse and its first two derivatives at one value, mirroring the
  R methods expression by expression,
  [`exp_floored()`](https://statmodels7.github.io/linkfunctions7/reference/exp_floored.md)’s
  derived guard included) and `lf7_clamp`
  ([`link_bounds_clamp()`](https://statmodels7.github.io/linkfunctions7/reference/link_bounds_clamp.md)
  on one value). A twin test holds each against the R method it stands
  for with [`identical()`](https://rdrr.io/r/base/identical.html).

## linkfunctions7 0.1.0

### Numerical behavior

- The derivatives of the transcendental links – logit, probit, cloglog,
  loglog, cauchit, rhobit and softplus – are computed by compiled
  kernels in `src/link_kernels.cpp`. The doubly bounded link and the
  softplus inverse reuse the logit kernel, scaled by the width and
  shifted one order. Measured through S7 dispatch,
  [`d4linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/d4linkinv.md)
  on the logit falls from 40 ms to 25 ms at a million points; what the
  port buys is not the transcendental, which costs the same in either
  language, but the vector temporaries an order-four polynomial
  allocates.

- The finite-difference weights and offsets come from `numericals7`,
  where the toolkit keeps one stencil library.

- [`linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/linkinv.md)
  returns a value strictly inside the open parameter interval. Nine of
  the eighteen links previously reached a bound or went non-finite under
  a large predictor: `linkinv(logit_link(), 37)` was exactly 1 and the
  round trip through
  [`linkfun()`](https://statmodels7.github.io/linkfunctions7/reference/linkfun.md)
  then `Inf`, which `distributions7` rejects, since it validates against
  open intervals. The correction is made in the generic, so every link
  inherits it and a user-written link needs no change. It fires only on
  exact equality with a bound; a value outside by a real margin is left
  alone, that being a complaint rather than a rounding artifact.

- [`cloglog_link()`](https://statmodels7.github.io/linkfunctions7/reference/cloglog_link.md)’s
  forward direction is computed through `log1p`. Written as
  `log(-log(1 - theta))` it returned `-Inf` for small `theta` where the
  value is finite and representable (-176.66 at `theta = 1.9e-77`),
  taking the four forward derivatives that divide by that logarithm with
  it.

- The floor the exponential links apply is derived from the quantity
  that binds rather than taken from a familiar constant. The old value,
  `.Machine$double.eps`, is the resolution of a double near one and not
  a lower bound on anything relevant; it corrupted `theta` for every
  `eta < -36`. The constraint is that the log link’s fourth derivative
  `-6/theta^4` must not overflow, giving `(24/double.xmax)^(1/4)`, about
  1.9e-77. The same correction applies to `cloglog` and the
  lower-bounded link.

- [`softplus_link()`](https://statmodels7.github.io/linkfunctions7/reference/softplus_link.md)
  is written in `u = -expm1(-a*theta)` rather than in `expm1(a*theta)`.
  The old form gave `NaN` for `d4linkfun` from `theta = 177` at `a = 1`
  and from `theta = 24` at `a = 30`, an ordinary value for a steep
  softplus.

### Derivatives

- A link that implements only
  [`linkfun()`](https://statmodels7.github.io/linkfunctions7/reference/linkfun.md)
  and
  [`linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/linkinv.md)
  is complete: the base class supplies the eight derivative generics
  numerically. Each fallback finds the highest order implemented
  analytically and applies one central stencil of order `k - m` to it,
  never a chain of first differences. On a log link defined with nothing
  the fourth derivative is accurate to about 2e-5; supplying the
  analytic first and second takes it to 2e-8.

- S7 classes are compared by name and package rather than by object
  identity. [`identical()`](https://rdrr.io/r/base/identical.html) on a
  class is `FALSE` for a class re-created from the same definition,
  which is what happens whenever the code is re-evaluated rather than
  loaded, so the base numerical fallback was taken for an analytic
  method and every fallback differentiated the order below it: the log
  link’s fourth derivative was wrong by a factor of 900.

### Validation and documentation

- [`check_link()`](https://statmodels7.github.io/linkfunctions7/reference/check_link.md)
  reports an unimplemented derivative as a failure. An order that raises
  left `NA` and broke the loop, and the summary reduced with
  `na.rm = TRUE`, so a link implementing only its first derivative was
  reported as passing all four orders.

- [`check_link()`](https://statmodels7.github.io/linkfunctions7/reference/check_link.md)
  distinguishes an order it verified from one it could not. An order
  supplied by a numerical fallback would be compared against a numerical
  differentiation of the order below, which is the same arithmetic twice
  and agrees however wrong the link is; such an order is left `NA` and
  printed as `[numerical]`.
  [`link_fallback_orders()`](https://statmodels7.github.io/linkfunctions7/reference/link_fallback_orders.md)
  answers the same question programmatically.

- Every object in the namespace has a help page, every exported topic a
  `\value` and an executable example, and `tests/testthat/test-docs.R`
  asks all of that plus whether the pkgdown index is complete.

- A vignette on working with link functions, a README with badges, a
  pkgdown site and continuous integration on five platforms.

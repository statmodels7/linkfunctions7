# linkfunctions7 0.5.2

* The inverse derivatives of `logit_link()`, and through the same kernel
  those of `bounded_link()` with two endpoints and of `softplus_link()`,
  are written in `p`, `1 - p` and their product, with `1 - p` computed from
  the predictor. Computed by subtraction from `p`, the complement lost its
  digits for a large positive predictor, and every order was 0 from
  `eta = 37`.

* The inverse derivatives of `rhobit_link()` are written in `tanh(eta)`
  and `sech(eta)^2`, with the latter computed from `exp(-2 * abs(eta))`.
  Computed as `1 - tanh(eta)^2`, every order lost its digits for a large
  `abs(eta)` and was 0 from `abs(eta) = 20`.

# linkfunctions7 0.5.1

* `plot()` on a link passes named graphical parameters to `graphics::plot()`
  for both panels, where they replace the defaults. They were documented as
  passed and were ignored. An unnamed argument in `...` signals an error.

* The documentation is rewritten in a plain register, and several stale
  statements are corrected: the constructors and `linkderiv()` and
  `linkinvderiv()` reach the fifth order, not the fourth; the class pages of
  the probit, log-log and power links count ten derivatives, not eight; the
  log-log link approaches 0 faster than 1; the softplus inverse approaches
  zero exponentially for a large negative predictor, not linearly; the
  square root link passes every check of `check_link()`; the pages of the
  second to fourth forward derivatives no longer describe the delta-method
  role of the first; and the forward fifth derivative of the exponential
  links is `Inf` below about `eta = -141`, not `-145`.

* The documentation of `linkinv()` and `link_bounds_clamp()` describes the
  clamp as implemented: only a value exactly on a bound, or an infinite
  value, is moved, and a value outside the bounds by more than rounding is
  returned unchanged (`linkinv(inverse_link(), -40)` is -0.025). Further
  corrections: the inverse derivative generics are computed from `eta`
  without the clamp; the forward derivative generics do not check the
  domain; `power_link(0)` returns the limit of the Box-Cox transformation,
  not of \eqn{\theta^\lambda}; every inverse derivative of an upper-bounded
  link changes sign; `eta_bounds()` is restricted for every non-zero power
  exponent; and a routed derivative call costs two to three times a direct
  one on a short vector.

# linkfunctions7 0.5.0

* `link_scalar_route()` returns the name and the parameters by which the
  compiled registry identifies a link, or `NULL`. Every link class has scalar
  C entries, `lf7_class_id()` and `lf7_inv12p()` beside `lf7_scalar_id()`,
  registered with `R_RegisterCCallable()`, which a compiled consumer such as
  the `gas()` filter of modelterms7 calls instead of calling back into R.
  The vector kernels and the entries share `src/link_points.h`.

# linkfunctions7 0.4.0

* Every link carries its fifth derivative analytically, in both directions:
  `d5linkfun()` and `d5linkinv()` are registered on all link classes,
  `linkderiv()` and `linkinvderiv()` route to order 5, and `check_link()`
  validates it. The order is needed by score-driven filters, where each order
  of differentiation of the predictor through its recursion involves one more
  order of the link and of the family: the curvature reaches the third, the
  directional third derivative the fourth, and the outer Hessian of a model
  with such a term the fifth.

  Each formula was derived from the order below it and checked against one
  Richardson pass on the analytic fourth derivative, with a negative control
  on the same grid (a coefficient 5 per cent out reads between 7.9e-03 and
  3.9, the shipped formula between 4.5e-12 and 7.4e-11). The higher orders
  follow from recurrences: the logistic polynomials satisfy
  \eqn{P_{k+1} = (1-2p)P_k + p(1-p)P_k'}, the probit's inverse derivatives are
  \eqn{(-1)^{k-1}He_{k-1}(\eta)\phi(\eta)}, the cloglog's satisfy
  \eqn{Q_{k+1} = (1-w)Q_k + wQ_k'}, and the softplus's forward numerators are
  the Eulerian numbers, 1, 11, 11, 1 at this order.

  In the inverse direction the softplus and the doubly bounded link reuse the
  logistic polynomials: the first one order down times a power of its
  steepness, the second scaled by the interval width. Tests assert both
  identities.

* A compiled kernel called with an unsupported order returns `NA`. Every
  kernel takes the order as an argument and previously reached its order-4
  branch through `default:`, so `lk_logit_inv_cpp(eta, 5L)` returned the
  fourth derivative. The routers have always rejected an unsupported order, so
  no public result was affected.

* The fifth forward derivative of the exponential links does not fit under
  the exponential floor, and the floor is unchanged. `exp_floor` is
  \eqn{(24/x_{\max})^{1/4}}, chosen so that \eqn{-6/\theta^4} cannot overflow;
  at the floor that derivative is -4.49e+307 against a `double.xmax` of
  1.80e+308. Covering the fifth would need \eqn{(120/x_{\max})^{1/5}} =
  5.82e-62, which would clamp \eqn{\theta} for \eqn{\eta} between about
  -176.7 and -141, where `linkinv()` is exact. `d5linkfun()` of a log, lower- or
  upper-bounded link is therefore `Inf` below about \eqn{\eta = -141}. In the
  lower tail the inverse direction, which the chain rule onto the link scale
  uses, is finite at every order, and a test covers both directions.

# linkfunctions7 0.3.0

* `logistic_deriv()` is tested against the compiled kernel of which it is
  the R statement. The logistic derivative polynomials are written twice, in
  R and as `logistic_poly()` in `src/link_kernels.cpp`, which the logit, the
  doubly bounded link and the softplus call. `test-logistic-twin.R` compares
  them at every order over 501 points and checks that each of the three links
  reaches the polynomial named in its description (unchanged, scaled by the
  interval width, and one order down times a power of the steepness).

  The comparison uses a tolerance, because both forms are Horner evaluations
  with multiply-adds that a compiler may fuse into FMA instructions. A
  negative control asserts that a polynomial wrong in one coefficient fails
  by more than 1e-3.

* `eta_bounds()` is exported. It returns the set of predictors that a link
  maps from, where `link_bounds` gives the set that it maps onto.
  `sqrt_link()`, `inverse_link()`, `inverse_sq_link()` and `power_link()` at
  any non-zero exponent map the positive half of the eta axis onto the positive
  half of the theta axis. A consumer that applies an inverse link to each
  coordinate of an unconstrained vector needs the map to be defined and
  injective on the whole real line, which fails outside those bounds:
  `linkinv(sqrt_link(), -2)` and `linkinv(sqrt_link(), 2)` are both 4.
  A package that depends on this one can now reject such a link at
  construction, where the error message can name the link.

# linkfunctions7 0.2.0

* Scalar C entry points for score-driven filters, registered with
  `R_RegisterCCallable`: `lf7_scalar_id` (identity and log; an unknown name
  returns -1 and the consumer keeps its R callbacks), `lf7_inv12` (the inverse
  and its first two derivatives at one value, matching the R methods
  expression by expression, including the guard of `exp_floored()`) and
  `lf7_clamp` (`link_bounds_clamp()` on one value). A test compares each with
  the corresponding R method using `identical()`.

# linkfunctions7 0.1.0

## Numerical behavior

* The derivatives of the transcendental links (logit, probit, cloglog,
  loglog, cauchit, rhobit and softplus) are computed by compiled kernels in
  `src/link_kernels.cpp`. The doubly bounded link and the softplus inverse
  reuse the logit kernel, scaled by the width and shifted one order. Measured
  through S7 dispatch, `d4linkinv()` on the logit takes 25 ms instead of
  40 ms at a million points; the saving comes from the vector temporaries
  that an order-four polynomial allocates, since the transcendental function
  costs the same in either language.

* The finite-difference weights and offsets come from `numericals7`.

* `linkinv()` returns a value strictly inside the open parameter interval.
  Previously several links reached a bound or returned a non-finite value for
  a large predictor: `linkinv(logit_link(), 37)` was exactly 1, so the round
  trip through `linkfun()` gave `Inf`, which `distributions7` rejects because
  it validates against open intervals. The clamp is applied in the generic,
  so every link inherits it, including a user-written one. It changes only a
  value exactly equal to a bound; a value outside by more than rounding is
  returned unchanged.

* `cloglog_link()`'s forward direction is computed through `log1p`. Written
  as `log(-log(1 - theta))` it returned `-Inf` for small `theta` where the
  value is finite (-176.66 at `theta = 1.9e-77`), and the four forward
  derivatives that divide by that logarithm with it.

* The floor of the exponential links is derived from the quantity that
  binds. The previous value, `.Machine$double.eps`, changed `theta` for every
  `eta < -36`. The constraint is that the log link's fourth derivative
  `-6/theta^4` must not overflow, which with a factor of four to spare gives
  `(24/double.xmax)^(1/4)`, about 1.9e-77. The same floor applies to `cloglog` and the lower-bounded link.

* `softplus_link()` is written in `u = -expm1(-a*theta)` instead of
  `expm1(a*theta)`. The previous form gave `NaN` for `d4linkfun` from
  `theta = 177` at `a = 1` and from `theta = 24` at `a = 30`.

## Derivatives

* A link that implements only `linkfun()` and `linkinv()` is complete: the
  base class supplies the derivative generics numerically. Each fallback
  takes the highest order implemented analytically and applies one central
  stencil of order `k - m` to it, never a chain of first differences. On a
  log link defined with nothing else, the fourth derivative is accurate to
  about 2e-5; supplying the analytic first and second derivatives brings it
  to 2e-8.

* S7 classes are compared by name and package, not by object identity.
  `identical()` on a class is `FALSE` for a class re-created from the same
  definition, which happens when the code is re-evaluated instead of loaded;
  the base numerical fallback was then taken for an analytic method, every
  fallback differentiated the order below it, and the log link's fourth
  derivative was wrong by a factor of 900.

## Validation and documentation

* `check_link()` reports an unimplemented derivative as a failure. An order
  that signaled an error was left `NA`, and the summary reduced with
  `na.rm = TRUE`, so a link implementing only its first derivative was
  reported as passing all four orders.

* `check_link()` distinguishes an order that it verified from one that it
  could not verify. An order supplied by a numerical fallback would be
  compared with a numerical differentiation of the order below, which is the
  same computation and agrees even for a wrong link; such an order is left
  `NA` and printed as `[numerical]`. `link_fallback_orders()` returns the
  same information.

* Every object in the namespace has a help page, every exported topic a
  `\value` and an executable example, and `tests/testthat/test-docs.R`
  checks all of this and whether the pkgdown index is complete.

* A vignette on working with link functions, a README with badges, a
  pkgdown site and continuous integration on five platforms.

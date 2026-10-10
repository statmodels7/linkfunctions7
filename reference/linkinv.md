# Evaluate Inverse Link Function

Maps a linear predictor back onto the parameter's domain.

## Usage

``` r
linkinv(x, eta)
```

## Arguments

- x:

  An object of class `link`.

- eta:

  A numeric vector of linear predictors.

## Value

A numeric vector of the same length as `eta`. A value that lands on a
bound or overflows is moved strictly inside `x@link_bounds`; a value
outside the bounds by more than rounding, `NA` and `NaN` are returned
unchanged.

## Details

A link is a bijection onto an open interval, so that a value returned by
the inverse link can be passed back to
[`linkfun()`](https://statmodels7.github.io/linkfunctions7/reference/linkfun.md),
or to a density that validates its parameters against open intervals. In
double precision the map is not quite onto: `plogis` is exactly 1 above
about \\\eta = 37\\, `lwr + exp(eta)` rounds to `lwr` once the
exponential falls below half an ulp of a non-zero `lwr`, both overflow
to infinity for large arguments, and the square-root, power, inverse and
inverse-square links reach a bound or overflow at \\\eta = 0\\.

The generic therefore clamps the result with
[`link_bounds_clamp()`](https://statmodels7.github.io/linkfunctions7/reference/link_bounds_clamp.md):
a value exactly equal to a finite bound is moved one or two units in the
last place inside it (to the smallest positive normal double when the
bound is zero), and an infinite value becomes the largest finite double
of that sign. A value outside the bounds by more than rounding, which
comes from a linear predictor outside
[`eta_bounds()`](https://statmodels7.github.io/linkfunctions7/reference/eta_bounds.md),
is returned unchanged, and so are `NA` and `NaN`. The clamp costs one
comparison per bound.

Because the clamp is applied in the generic, every link inherits it,
including a user-defined one, and a method can be written as the plain
formula without a guard of its own.

## See also

[`linkfun()`](https://statmodels7.github.io/linkfunctions7/reference/linkfun.md),
[`link_bounds_clamp()`](https://statmodels7.github.io/linkfunctions7/reference/link_bounds_clamp.md)

## Examples

``` r
linkinv(logit_link(), c(-1, 0, 1))
#> [1] 0.2689414 0.5000000 0.7310586
linkinv(log_link(), c(0, 1))
#> [1] 1.000000 2.718282

# far out, where the arithmetic saturates: strictly inside (0, 1) rather
# than exactly 1, so the round trip still returns a number
linkinv(logit_link(), 40) < 1
#> [1] TRUE
linkfun(logit_link(), linkinv(logit_link(), 40))
#> [1] 36.04365
```

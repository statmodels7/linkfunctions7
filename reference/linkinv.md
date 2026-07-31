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

A numeric vector, strictly inside `x@link_bounds`.

## Details

**The result is always strictly inside `link_bounds`**, and that is
enforced here rather than left to each method. A link is a bijection
onto an *open* interval, which is exactly what makes it useful — a value
it returns can be handed to
[`linkfun`](https://statmodels7.github.io/linkfunctions7/reference/linkfun.md)
and come back, or to a density that validates its parameters against
open intervals. In double precision the bijection is not quite onto:
`plogis` is exactly 1 above about \\\eta = 37\\, `lwr + exp(eta)` rounds
to `lwr` once the exponential falls below half an ulp of a non-zero
`lwr`, and both overflow to infinity far out. Nine of the links shipped
here reach a bound somewhere in \\\lvert \eta \rvert \le 800\\.

So the generic clamps: a result at or beyond a finite bound becomes the
nearest double strictly inside it, and a non-finite result becomes the
largest finite double of that sign. Both are derived rather than chosen
— they are the extremes of what "strictly inside, and a usable number"
permits, and no tolerance is invented. The clamp costs one comparison
per bound and fires only in the tails.

Doing it in the generic body means every link inherits it, including one
you write yourself, and that a method can be written as the plain
mathematical formula without a guard of its own.

## See also

[`linkfun`](https://statmodels7.github.io/linkfunctions7/reference/linkfun.md),
[`link_bounds_clamp`](https://statmodels7.github.io/linkfunctions7/reference/link_bounds_clamp.md)

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

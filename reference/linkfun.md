# Evaluate Forward Link Function

Carries a parameter from its own domain onto the real line, \\\eta =
g(\theta)\\.

## Usage

``` r
linkfun(x, theta)
```

## Arguments

- x:

  An object of class `link`.

- theta:

  A numeric vector of parameters.

## Value

A numeric vector of the linear predictor.

## Details

A link is a strictly monotone differentiable bijection \\g : \Theta \to
\mathbb{R}\\ from the parameter's open domain `x@link_bounds` onto the
whole line, so that an unconstrained optimizer may work in \\\eta\\
while \\\theta = g^{-1}(\eta)\\ stays admissible at every point.
[`linkinv`](https://statmodels7.github.io/linkfunctions7/reference/linkinv.md)
evaluates that inverse and is the only other method a link has to
supply; the eight derivative generics have numerical fallbacks derived
from the pair.

## Examples

``` r
linkfun(logit_link(), c(0.25, 0.5, 0.75))
#> [1] -1.098612  0.000000  1.098612
linkfun(log_link(), c(1, exp(1)))
#> [1] 0 1
```

# Clamp a Parameter Strictly Inside Its Domain

Moves a value that has reached or passed a bound to the nearest double
strictly inside it, and a non-finite value to the largest finite double
of that sign. Applied by
[`linkinv`](https://statmodels7.github.io/linkfunctions7/reference/linkinv.md)
to every link.

## Usage

``` r
link_bounds_clamp(theta, bounds)
```

## Arguments

- theta:

  A numeric vector, as a method computed it.

- bounds:

  The link's `link_bounds`, a length-2 numeric vector.

## Value

`theta`, with any value that has landed exactly on a bound moved just
inside it and any infinity brought back to the largest finite double.
`NA` and `NaN` pass through untouched, as does a value strictly outside
by a real margin: converting either would hide something rather than fix
it.

## Details

A link is documented as a bijection onto an **open** interval, and in
exact arithmetic it is. In double precision it is not: `plogis(37)` is
exactly 1, `2 + exp(-40)` is exactly 2, and `exp(800)` is infinite. A
caller then receives a probability of exactly 1, or a variance of
exactly 0, and the next thing they do is take its logarithm or divide by
it.

The correction is the smallest one that can work and it is derived, not
chosen. What binds is "strictly inside, and a number you can compute
with", and its two extremes are the neighbouring representable double
and the largest finite one.

### Why the bump is a relative one

R has no `nextafter`, and the arithmetic substitute has to respect that
**the spacing of doubles is absolute near a non-zero bound**. One ulp at
2 is about 4.4e-16 while one ulp at 1e-300 is about 1e-316, so a single
additive constant cannot serve both. `b + |b| * eps` is one to two ulps
from `b` at any magnitude, which is strictly inside and as close as
arithmetic reliably gets.

A bound at zero is the exception and needs no bump, since there the
spacing is relative all the way down to 1e-308 and the exponential links
already floor at
[`exp_floor`](https://statmodels7.github.io/linkfunctions7/reference/exp_floor.md).
The clamp therefore leaves an exact zero bound to
[`exp_floor`](https://statmodels7.github.io/linkfunctions7/reference/exp_floor.md)
and uses the smallest positive normal only if something has still landed
on it.

## See also

[`linkinv`](https://statmodels7.github.io/linkfunctions7/reference/linkinv.md)

## Examples

``` r
link_bounds_clamp(c(0, 0.5, 1), c(0, 1))
#> [1] 2.225074e-308  5.000000e-01  1.000000e+00
link_bounds_clamp(c(2, 3, Inf), c(2, Inf))
#> [1]  2.000000e+00  3.000000e+00 1.797693e+308
```

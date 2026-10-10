# Clamp a Parameter Strictly Inside Its Domain

Moves a value that lies exactly on a finite bound to a double strictly
inside it, one or two units in the last place away, and an infinite
value to the largest finite double of that sign. Applied by
[`linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/linkinv.md)
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
`NA` and `NaN` pass through unchanged, as does a value outside the
bounds by more than rounding, since changing either would hide an error
in the calling code.

## Details

In exact arithmetic a link is a bijection onto an open interval, but in
double precision the endpoints can be reached: `plogis(37)` is exactly
1, `2 + exp(-40)` is exactly 2, and `exp(800)` is infinite. A caller
would then receive a probability of exactly 1 or a variance of exactly
0, and taking its logarithm or dividing by it fails.

The clamp moves such a value by the smallest amount that makes it
strictly inside the interval and finite: to a neighboring representable
double, or to the largest finite double.

### The size of the step

R has no `nextafter()`, so the step is computed arithmetically. Near a
non-zero bound \\b\\ the spacing of doubles is set by the magnitude of
\\b\\: one ulp at 2 is about 4.4e-16, while one ulp at 1e-300 is about
1e-316, so no single additive constant serves both. The value
`b + |b| * eps` lies one to two ulps from `b` at any magnitude, which is
strictly inside the interval and as close to `b` as the arithmetic
reliably allows.

A bound at zero needs no such step, because the exponential links
already floor their result at
[`exp_floor()`](https://statmodels7.github.io/linkfunctions7/reference/exp_floor.md).
A value that still lands exactly on a zero bound is moved to the
smallest positive normal double.

## See also

[`linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/linkinv.md)

## Examples

``` r
link_bounds_clamp(c(0, 0.5, 1), c(0, 1))
#> [1] 2.225074e-308  5.000000e-01  1.000000e+00
link_bounds_clamp(c(2, 3, Inf), c(2, Inf))
#> [1]  2.000000e+00  3.000000e+00 1.797693e+308
```

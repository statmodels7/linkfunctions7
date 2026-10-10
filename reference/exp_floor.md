# The Floor of the Exponential Links

The floor applied to `exp(eta)` by every link whose inverse is an
exponential
([`log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.md),
[`cloglog_link()`](https://statmodels7.github.io/linkfunctions7/reference/cloglog_link.md),
and the lower- and upper-bounded links).

## Usage

``` r
exp_floor
```

## Format

A length-one numeric vector.

## Value

A length-one numeric vector, about `1.9e-77`, at which \\-6/\theta^4\\
evaluates to `-4.5e307`.

## Details

The floor exists so that a parameter returned as \\\theta\\ can be
divided into without producing `Inf`: the forward derivatives of these
links are \\1/\theta\\, \\-1/\theta^2\\, \\2/\theta^3\\ and
\\-6/\theta^4\\, and the fourth is the binding one. Solving \\6/\theta^4
\le\\ `double.xmax` and keeping a factor of four in hand gives
`(24 / .Machine$double.xmax)^0.25`, about `1.9e-77`, at which
\\-6/\theta^4\\ evaluates to `-4.5e307`.

The floor is the lowest value that this constraint allows, and it keeps
\\\theta\\ exact down to \\\eta \approx -176.7\\. The fifth forward
derivative, \\24/\theta^5\\, does not fit under it and is `Inf` for
\\\eta\\ below about \\-141\\.

## See also

[`exp_floored()`](https://statmodels7.github.io/linkfunctions7/reference/exp_floored.md)

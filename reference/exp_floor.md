# The Smallest Parameter Value the Exponential Links Will Report

The floor applied to `exp(eta)` by every link whose inverse is an
exponential
([`log_link`](https://statmodels7.github.io/linkfunctions7/reference/log_link.md),
[`cloglog_link`](https://statmodels7.github.io/linkfunctions7/reference/cloglog_link.md),
and the lower- and upper-bounded links).

## Usage

``` r
exp_floor
```

## Format

A length-one numeric vector.

## Details

The floor exists so that a parameter reported as \\\theta\\ can be
divided into without producing `Inf`: the forward derivatives of these
links are \\1/\theta\\, \\-1/\theta^2\\, \\2/\theta^3\\ and
\\-6/\theta^4\\, and the fourth is the binding one. Solving \\6/\theta^4
\le\\ `double.xmax` and keeping a factor of four in hand gives
`(24 / .Machine$double.xmax)^0.25`, about `1.9e-77`, at which
\\-6/\theta^4\\ evaluates to `-4.5e307`.

The point of choosing it this way is that the floor should be as *low*
as that constraint allows, not as high as seems safe. It was previously
`.Machine$double.eps`, which is 61 orders of magnitude higher than
necessary and silently corrupted \\\theta\\ for every \\\eta \< -36\\:
`linkinv(log_link(), -40)` returned `2.2e-16` instead of `4.2e-18`, and
the round trip came back `-36.04` instead of `-40`. The present value
keeps \\\theta\\ exact down to \\\eta \approx -177\\ while leaving every
derivative just as finite as before.

## See also

[`exp_floored`](https://statmodels7.github.io/linkfunctions7/reference/exp_floored.md)

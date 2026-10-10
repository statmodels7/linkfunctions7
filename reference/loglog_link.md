# The Log-Log Link Function

The log-log link \\\eta = -\log(-\log\theta)\\ on \\(0, 1)\\, with
inverse \\\theta = \exp(-e^{-\eta})\\; the mirror image of
[`cloglog_link()`](https://statmodels7.github.io/linkfunctions7/reference/cloglog_link.md).

## Usage

``` r
loglog_link()
```

## Value

An S7 object of class `LogLogLink`, inheriting from
[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
whose methods compute the link, its inverse and their derivatives to the
fifth order.

## Details

The log-log link is \\\eta = -\log(-\log\theta)\\, with inverse \\\theta
= \exp(-\exp(-\eta))\\.

Unlike the logit and the probit, the link is asymmetric: \\\theta\\
approaches 0 faster than 1, which is the mirror image of
[`cloglog_link()`](https://statmodels7.github.io/linkfunctions7/reference/cloglog_link.md).
The domain of \\\theta\\ is \\(0, 1)\\.

## See also

[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
[`cloglog_link()`](https://statmodels7.github.io/linkfunctions7/reference/cloglog_link.md),
[`logit_link()`](https://statmodels7.github.io/linkfunctions7/reference/logit_link.md)

## Examples

``` r
lk <- loglog_link()
lk
#> S7 Link Object: loglog
#>   - Parameter domain (theta): (0, 1)

p <- c(0.1, 0.5, 0.9)
eta <- linkfun(lk, p)
eta
#> [1] -0.8340324  0.3665129  2.2503673
linkinv(lk, eta)
#> [1] 0.1 0.5 0.9

# the mirror image of the cloglog link
linkinv(loglog_link(), 1)
#> [1] 0.6922006
1 - linkinv(cloglog_link(), -1)
#> [1] 0.6922006

linkderiv(lk, 0.5, order = 3)
#> [1] 21.17476
```

# The Log-Log Link Function

Creates an S7 object of class `link` implementing the Log-Log
transformation. This link is the asymmetric opposite of the `cloglog`
link and is often used in survival analysis (e.g., Gompertz models) or
for modeling binary data with specific asymmetry requirements.

## Usage

``` r
loglog_link()
```

## Value

An S7 object of class `LogLogLink` (inheriting from `link`) containing
the transformation functions and their exact analytical derivatives up
to the fourth order.

## Details

The Log-Log link is mathematically defined as \\\eta =
-\log(-\log(\theta))\\. Consequently, the inverse link is derived as
\\\theta = \exp(-\exp(-\eta))\\.

**Asymmetry:** Unlike the symmetric Logit or Probit links, the Log-Log
link is highly asymmetric. It is particularly suitable for modeling
events where the probability approaches 0 very slowly but approaches 1
sharply. This represents the mathematical reverse of the `cloglog` link.

The strict mathematical domain for \\\theta\\ is `c(0, 1)`.

## See also

[`link`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
[`cloglog_link`](https://statmodels7.github.io/linkfunctions7/reference/cloglog_link.md),
[`logit_link`](https://statmodels7.github.io/linkfunctions7/reference/logit_link.md)

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

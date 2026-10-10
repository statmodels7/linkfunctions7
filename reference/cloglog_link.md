# The Complementary Log-Log (ClogLog) Link Function

The complementary log-log link \\\eta = \log(-\log(1-\theta))\\ on \\(0,
1)\\, with inverse \\\theta = 1 - \exp(-e^\eta)\\; asymmetric about
\\\theta = 1/2\\.

## Usage

``` r
cloglog_link()
```

## Value

An S7 object of class `ClogLogLink`, inheriting from
[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
whose methods compute the link, its inverse and their derivatives to the
fifth order.

## Details

The complementary log-log link is \\\eta = \log(-\log(1 - \theta))\\,
with inverse \\\theta = 1 - \exp(-\exp(\eta))\\.

Unlike the logit and the probit, the link is asymmetric: \\\theta\\
approaches 1 faster than 0. The inverse link is the distribution
function of the Gumbel distribution for minima. The link is used in
discrete-time survival analysis, where it gives a proportional-hazards
model, and for rare events.

The domain of \\\theta\\ is \\(0, 1)\\.

## See also

[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
[`logit_link()`](https://statmodels7.github.io/linkfunctions7/reference/logit_link.md),
[`loglog_link()`](https://statmodels7.github.io/linkfunctions7/reference/loglog_link.md)

## Examples

``` r
lk <- cloglog_link()
lk
#> S7 Link Object: cloglog
#>   - Parameter domain (theta): (0, 1)

p <- c(0.1, 0.5, 0.9)
eta <- linkfun(lk, p)
eta
#> [1] -2.2503673 -0.3665129  0.8340324
linkinv(lk, eta)
#> [1] 0.1 0.5 0.9

# asymmetric: it reaches 1 sharply and 0 slowly, the mirror of loglog
linkinv(cloglog_link(), c(-2, 2))
#> [1] 0.126577 0.999382
linkinv(loglog_link(),  c(-2, 2))
#> [1] 0.000617979 0.873423018

d2linkinv(lk, 0)
#> [1] 0
```

# The Cauchit Link Function

The cauchit link \\\eta = \tan(\pi(\theta - 1/2))\\ on \\(0, 1)\\, the
Cauchy quantile function; heavier-tailed than the logit or the probit.

## Usage

``` r
cauchit_link()
```

## Value

An S7 object of class `CauchitLink`, inheriting from
[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
whose methods compute the link, its inverse and their derivatives to the
fifth order.

## Details

The cauchit link is \\\eta = \tan(\pi(\theta - 1/2))\\, computed as
`qcauchy(theta)`. The inverse link is the standard Cauchy distribution
function \\\theta = \arctan(\eta)/\pi + 1/2\\, computed as
`pcauchy(eta)`.

The Cauchy distribution has heavier tails than the logistic and the
normal, so the probability approaches 0 and 1 more slowly than under the
logit or the probit link. The link is therefore less sensitive to
observations whose linear predictor is extreme.

The domain of \\\theta\\ is \\(0, 1)\\.

## See also

[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
[`logit_link()`](https://statmodels7.github.io/linkfunctions7/reference/logit_link.md),
[`probit_link()`](https://statmodels7.github.io/linkfunctions7/reference/probit_link.md)

## Examples

``` r
lk <- cauchit_link()
lk
#> S7 Link Object: cauchit
#>   - Parameter domain (theta): (0, 1)

p <- c(0.1, 0.5, 0.9)
eta <- linkfun(lk, p)
eta
#> [1] -3.077684  0.000000  3.077684
linkinv(lk, eta)
#> [1] 0.1 0.5 0.9

# heavy tails: the same eta is far less extreme than under a logit
linkinv(cauchit_link(), 5)
#> [1] 0.937167
linkinv(logit_link(), 5)
#> [1] 0.9933071

dlinkinv(lk, 0)            # 1 / pi
#> [1] 0.3183099
```

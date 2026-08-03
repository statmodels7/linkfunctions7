# The Cauchit Link Function

The cauchit link \\\eta = \tan(\pi(\theta - 1/2))\\ on \\(0, 1)\\, the
Cauchy quantile function; heavier-tailed than the logit or the probit.

## Usage

``` r
cauchit_link()
```

## Value

An S7 object of class `CauchitLink` (inheriting from `link`) containing
the transformation functions and their exact analytical derivatives up
to the fourth order.

## Details

The Cauchit link is defined mathematically as \\\eta = \tan(\pi(\theta -
0.5))\\, which corresponds perfectly to `qcauchy(theta)`. The inverse
link is the standard Cauchy CDF \\\theta = \frac{1}{\pi} \arctan(\eta) +
0.5\\, computed via `pcauchy(eta)`.

**Heavy Tails:** Unlike the Logit or Probit links, the Cauchit link has
exceedingly heavier tails. This makes it particularly robust and useful
for modeling binary data where the probability approaches 0 or 1 very
slowly, or when the dataset contains severe outliers that might
disproportionately influence the fit of light-tailed link functions.

The strictly valid mathematical domain for \\\theta\\ is `c(0, 1)`.

## See also

[`link`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
[`logit_link`](https://statmodels7.github.io/linkfunctions7/reference/logit_link.md),
[`probit_link`](https://statmodels7.github.io/linkfunctions7/reference/probit_link.md)

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

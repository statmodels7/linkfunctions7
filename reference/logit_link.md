# The Logit Link Function

The logit link \\\eta = \log(\theta/(1-\theta))\\ on \\(0, 1)\\, with
inverse \\\theta = 1/(1+e^{-\eta})\\; the canonical link for a
probability, whose linear predictor is the log-odds.

## Usage

``` r
logit_link()
```

## Value

An S7 object of class `LogitLink` (inheriting from `link`) containing
the transformation functions and their exact analytical derivatives up
to the fourth order.

## Details

The Logit link is defined mathematically as \\\eta =
\log(\frac{\theta}{1 - \theta})\\. The inverse link is the standard
logistic function (sigmoid): \\\theta = \frac{1}{1 + \exp(-\eta)}\\.

The link is symmetric about \\\theta = 0.5\\, where \\\eta = 0\\, and
the linear predictor is the log-odds of the event probability. The
domain of \\\theta\\ is \\(0, 1)\\.

The implementation delegates to
[`stats::qlogis`](https://rdrr.io/r/stats/Logistic.html) and
[`stats::plogis`](https://rdrr.io/r/stats/Logistic.html), which remain
accurate near both boundaries.

## See also

[`link`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
[`probit_link`](https://statmodels7.github.io/linkfunctions7/reference/probit_link.md),
[`cloglog_link`](https://statmodels7.github.io/linkfunctions7/reference/cloglog_link.md)

## Examples

``` r
lk <- logit_link()
lk
#> S7 Link Object: logit
#>   - Parameter domain (theta): (0, 1)

p <- c(0.1, 0.5, 0.9)
eta <- linkfun(lk, p)      # log-odds
eta
#> [1] -2.197225  0.000000  2.197225
linkinv(lk, eta)
#> [1] 0.1 0.5 0.9

# the inverse derivatives are polynomials in p: the first is the variance
# of a Bernoulli, p(1 - p)
dlinkinv(lk, 0)
#> [1] 0.25

# all four orders at once
vapply(1:4, function(k) linkinvderiv(lk, 0, order = k), numeric(1))
#> [1]  0.250  0.000 -0.125  0.000

# every mathematical property is checkable
check_link(lk)
#> Checking S7 Link Object: logit 
#>   [1] Invertibility (Theta space): [PASSED] 
#>   [2] Invertibility (Eta space):   [PASSED] 
#>   [3] Strict Monotonicity:         [PASSED] 
#>   [4] Inverse Function Theorem:    [PASSED] 
#>   [5] Link Derivatives:            [PASSED] 
#>   [6] Inverse Link Derivatives:    [PASSED] 
```

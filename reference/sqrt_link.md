# The Square Root Link Function

The square-root link \\\eta = \sqrt{\theta}\\ on \\(0, \infty)\\; its
image is \\(0, \infty)\\.

## Usage

``` r
sqrt_link()
```

## Value

An S7 object of class `SqrtLink`, inheriting from
[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
whose methods compute the link, its inverse and their derivatives to the
fifth order.

## Details

The square-root link is \\\eta = \sqrt{\theta}\\, with inverse \\\theta
= \eta^2\\.

The inverse \\\eta^2\\ is defined for a negative \\\eta\\ as well, but
there it is not one-to-one, since \\\eta\\ and \\-\eta\\ give the same
\\\theta\\. The linear predictor is therefore meant to stay positive;
[`eta_bounds()`](https://statmodels7.github.io/linkfunctions7/reference/eta_bounds.md)
returns this range.

The domain of \\\theta\\ is \\(0, \infty)\\.

## See also

[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
[`power_link()`](https://statmodels7.github.io/linkfunctions7/reference/power_link.md),
[`log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.md)

## Examples

``` r
lk <- sqrt_link()
lk
#> S7 Link Object: sqrt
#>   - Parameter domain (theta): (0, Inf)

theta <- c(0.25, 1, 4)
eta <- linkfun(lk, theta)
eta
#> [1] 0.5 1.0 2.0
linkinv(lk, eta)
#> [1] 0.25 1.00 4.00

# the inverse is a quadratic, so the third and fourth derivatives vanish
d2linkinv(lk, c(1, 2))
#> [1] 2 2
d3linkinv(lk, c(1, 2))
#> [1] 0 0

# the same link as the power family at lambda = 1/2
linkfun(power_link(0.5), 4)
#> [1] 2
```

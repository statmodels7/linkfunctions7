# The Inverse (Reciprocal) Link Function

The reciprocal link \\\eta = 1/\theta\\ on \\(0, \infty)\\, the
canonical link of the Gamma family; its image is \\(0, \infty)\\, not
the whole real line.

## Usage

``` r
inverse_link()
```

## Value

An S7 object of class `InverseLink`, inheriting from
[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
whose methods compute the link, its inverse and their derivatives to the
fifth order.

## Details

The inverse link is \\\eta = 1/\theta\\, and it is its own inverse,
\\\theta = 1/\eta\\.

It is the canonical link of the Gamma family, used when the mean is
modeled as the reciprocal of a linear predictor.

The domain of \\\theta\\ is \\(0, \infty)\\, and so is the image of the
link, so the linear predictor must stay positive during optimization;
[`eta_bounds()`](https://statmodels7.github.io/linkfunctions7/reference/eta_bounds.md)
returns this range.

## See also

[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
[`identity_link()`](https://statmodels7.github.io/linkfunctions7/reference/identity_link.md)

## Examples

``` r
lk <- inverse_link()
lk
#> S7 Link Object: inverse
#>   - Parameter domain (theta): (0, Inf)

theta <- c(0.5, 1, 2)
eta <- linkfun(lk, theta)
eta
#> [1] 2.0 1.0 0.5
linkinv(lk, eta)           # the map is its own inverse
#> [1] 0.5 1.0 2.0

dlinkfun(lk, theta)
#> [1] -4.00 -1.00 -0.25

# the canonical link for a Gamma mean; note eta must keep one sign
linkinv(lk, c(0.5, 2))
#> [1] 2.0 0.5
```

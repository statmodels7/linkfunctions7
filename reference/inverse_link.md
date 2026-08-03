# The Inverse (Reciprocal) Link Function

The reciprocal link \\\eta = 1/\theta\\ on \\(0, \infty)\\, the
canonical link of the Gamma family; its image is \\(0, \infty)\\, not
the whole real line.

## Usage

``` r
inverse_link()
```

## Value

An S7 object of class `InverseLink` (inheriting from `link`) containing
the transformation functions and their exact analytical derivatives up
to the fourth order.

## Details

The Inverse link is defined as \\\eta = 1/\theta\\. The inverse link
function is therefore perfectly symmetric: \\\theta = 1/\eta\\.

This link is typically used for modeling positive continuous data where
the mean is inversely proportional to the linear predictor (e.g., in
Gamma regression).

The domain of \\\theta\\ is conventionally `c(0, Inf)`. Care must be
taken to ensure the linear predictor \\\eta\\ remains strictly positive
(or strictly negative) during optimization to avoid division by zero or
mapping to invalid negative parameter values.

## See also

[`link`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
[`identity_link`](https://statmodels7.github.io/linkfunctions7/reference/identity_link.md)

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

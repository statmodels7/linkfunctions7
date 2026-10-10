# The Inverse Square Link Function

The inverse-square link \\\eta = 1/\theta^2\\ on \\(0, \infty)\\, the
canonical link of the inverse Gaussian family; its image is \\(0,
\infty)\\.

## Usage

``` r
inverse_sq_link()
```

## Value

An S7 object of class `InverseSqLink`, inheriting from
[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
whose methods compute the link, its inverse and their derivatives to the
fifth order.

## Details

The inverse-square link is \\\eta = 1/\theta^2\\, with inverse \\\theta
= 1/\sqrt{\eta}\\.

It is the canonical link of the inverse Gaussian family, whose variance
is proportional to the cube of the mean.

The domain of \\\theta\\ is \\(0, \infty)\\, and so is the image of the
link. The inverse link and its derivatives return `NaN` for a negative
\\\eta\\, so the linear predictor must stay positive during
optimization;
[`eta_bounds()`](https://statmodels7.github.io/linkfunctions7/reference/eta_bounds.md)
returns this range.

## See also

[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
[`inverse_link()`](https://statmodels7.github.io/linkfunctions7/reference/inverse_link.md)

## Examples

``` r
lk <- inverse_sq_link()
lk
#> S7 Link Object: inverse_sq
#>   - Parameter domain (theta): (0, Inf)

theta <- c(0.5, 1, 2)
eta <- linkfun(lk, theta)  # 1 / theta^2
eta
#> [1] 4.00 1.00 0.25
linkinv(lk, eta)
#> [1] 0.5 1.0 2.0

# the canonical link of the inverse Gaussian; eta must stay positive
dlinkinv(lk, c(0.5, 1, 4))
#> [1] -1.414214 -0.500000 -0.062500
```

# The Logarithmic Link Function

The log link \\\eta = \log\theta\\ on \\(0, \infty)\\, with inverse
\\\theta = e^\eta\\; the canonical link for a positive parameter.

## Usage

``` r
log_link()
```

## Value

An S7 object of class `LogLink` (inheriting from `link`) containing the
transformation functions and their exact analytical derivatives up to
the fourth order.

## Details

The Log link is defined mathematically as \\\eta = \log(\theta)\\. The
inverse link is the exponential function \\\theta = \exp(\eta)\\.

The inverse function of this link is its own derivative. Therefore, the
parameter \\\theta\\ and all its derivatives with respect to \\\eta\\
are equal to \\\exp(\eta)\\.

The valid mathematical domain of \\\theta\\ is `c(0, Inf)`.

**Numerical Stability:** The inverse link and its derivatives are
bounded below by `exp_floor`, which is `.Machine$double.xmin^0.25`,
about `1.2e-77`. This prevents underflow to exactly zero for large
negative \\\eta\\, which would produce `Inf` when the forward
derivatives divide by \\\theta\\; the fourth of them divides by
\\\theta^4\\, and that is what sets the value. The floor is low enough
that \\\theta\\ is exact down to \\\eta \approx -177\\.

## See also

[`link`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
[`inverse_link`](https://statmodels7.github.io/linkfunctions7/reference/inverse_link.md)

## Examples

``` r
lk <- log_link()
lk
#> S7 Link Object: log
#>   - Parameter domain (theta): (0, Inf)

theta <- c(0.5, 1, 10)
eta <- linkfun(lk, theta)
eta
#> [1] -0.6931472  0.0000000  2.3025851
linkinv(lk, eta)
#> [1]  0.5  1.0 10.0

# the exponential is its own derivative, so every inverse derivative agrees
dlinkinv(lk, eta)
#> [1]  0.5  1.0 10.0
d4linkinv(lk, eta)
#> [1]  0.5  1.0 10.0

# forward derivatives to fourth order
linkderiv(lk, theta, order = 4)
#> [1] -96.0000  -6.0000  -0.0006
```

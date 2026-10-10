# The Logarithmic Link Function

The log link \\\eta = \log\theta\\ on \\(0, \infty)\\, with inverse
\\\theta = e^\eta\\; the canonical link for a positive parameter.

## Usage

``` r
log_link()
```

## Value

An S7 object of class `LogLink`, inheriting from
[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
whose methods compute the link, its inverse and their derivatives to the
fifth order.

## Details

The log link is \\\eta = \log(\theta)\\, and its inverse is the
exponential function \\\theta = \exp(\eta)\\.

The exponential is its own derivative, so every derivative of the
inverse link with respect to \\\eta\\ equals \\\exp(\eta)\\.

The domain of \\\theta\\ is \\(0, \infty)\\.

The inverse link and its derivatives are bounded below by
[`exp_floor()`](https://statmodels7.github.io/linkfunctions7/reference/exp_floor.md),
\\(24/x\_{\max})^{1/4} \approx 1.9 \times 10^{-77}\\, where
\\x\_{\max}\\ is the largest finite double. The floor prevents an
underflow to exactly zero for a large negative \\\eta\\, which would
make the forward derivatives infinite; its value is set by the fourth of
them, which divides by \\\theta^4\\. The floor is low enough that
\\\theta\\ is exact down to \\\eta \approx -176.7\\.

## See also

[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
[`inverse_link()`](https://statmodels7.github.io/linkfunctions7/reference/inverse_link.md)

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

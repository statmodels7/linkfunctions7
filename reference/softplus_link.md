# The Softplus Link Function

The softplus link with scale \\a\\: the inverse is \\\theta = \log(1 +
e^{a\eta})/a\\, a smooth approximation of \\\max(0, \eta)\\ that
sharpens as \\a\\ grows.

## Usage

``` r
softplus_link(a = 1)
```

## Arguments

- a:

  The scale parameter, which sets the steepness: a single number
  strictly greater than 0. Defaults to 1.

## Value

An S7 object of class `SoftplusLink`, inheriting from
[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
whose methods compute the link, its inverse and their derivatives to the
fifth order; its `link_params` holds `a`.

## Details

The softplus link expresses the parameter \\\theta\\ as the softplus of
the linear predictor \\\eta\\:

- inverse link: \\\theta = \frac{1}{a} \log(1 + \exp(a \eta))\\;

- link: \\\eta = \frac{1}{a} \log(\exp(a \theta) - 1)\\.

For a large negative \\\eta\\, \\\theta \approx 0\\. For a large
positive \\\eta\\, \\\theta \approx \eta\\, which grows linearly where
the inverse of a log link would grow as \\\exp(\eta)\\.

Both directions are written so that no intermediate quantity grows with
\\a\theta\\ or \\a\eta\\. The inverse link uses the log-sum-exp form,
and the forward link and its derivatives are expressed in \\u = 1 -
e^{-a\theta}\\ instead of \\e^{a\theta} - 1\\. The latter overflows once
\\a\theta\\ exceeds about 709, and earlier at the higher orders, because
the derivatives divide by a power of it.

The domain of \\\theta\\ is \\(0, \infty)\\.

## See also

[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
[`log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.md),
[`identity_link()`](https://statmodels7.github.io/linkfunctions7/reference/identity_link.md)

## Examples

``` r
lk <- softplus_link(a = 2)
lk
#> S7 Link Object: softplus(a=2)
#>   - Parameter domain (theta): (0, Inf)
#>   - Link parameters: a = 2

theta <- c(0.5, 1, 5)
eta <- linkfun(lk, theta)
eta
#> [1] 0.2706624 0.9272933 4.9999773
linkinv(lk, eta)          # back to theta
#> [1] 0.5 1.0 5.0

# derivatives of either direction, to fourth order
dlinkfun(lk, theta)
#> [1] 1.581977 1.156518 1.000045
d4linkinv(lk, eta)
#> [1] -0.735332435  0.278864491  0.000363084

# unlike the log link, softplus is asymptotically linear in eta
linkinv(softplus_link(), c(1, 10, 100))
#> [1]   1.313262  10.000045 100.000000
```

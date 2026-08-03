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

  A numeric value specifying the scaling parameter
  (smoothness/steepness). Must be strictly positive. Defaults to 1.

## Value

An S7 object of class `SoftplusLink` (inheriting from `link`) containing
the transformation functions, their exact analytical derivatives up to
the fourth order, and the parameter `a`.

## Details

The Softplus link describes the relationship where the response
parameter \\\theta\\ is the Softplus of the linear predictor \\\eta\\.

Mathematically:

- Inverse Link (Softplus): \\\theta = \frac{1}{a} \log(1 + \exp(a
  \eta))\\

- Link Function: \\\eta = \frac{1}{a} \log(\exp(a \theta) - 1)\\

**Behavior:** For large negative \\\eta\\, \\\theta \approx 0\\. For
large positive \\\eta\\, \\\theta \approx \eta\\ (linear behavior),
whereas a Log link would imply \\\theta = \exp(\eta)\\ (exponential
behavior).

**Numerical Stability:** Both directions are written so that no
intermediate quantity grows with \\a\theta\\ or \\a\eta\\. The inverse
link uses the log-sum-exp form, and the forward link and its derivatives
are expressed in \\u = 1 - e^{-a\theta}\\ rather than in \\e^{a\theta} -
1\\, which overflows once \\a\theta\\ passes about 709 — and, because
the derivatives divide by its fourth power, well before that at the
higher orders.

The mathematical domain of \\\theta\\ is `c(0, Inf)`.

## See also

[`link`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
[`log_link`](https://statmodels7.github.io/linkfunctions7/reference/log_link.md),
[`identity_link`](https://statmodels7.github.io/linkfunctions7/reference/identity_link.md)

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

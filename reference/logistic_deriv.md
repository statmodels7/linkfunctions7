# Derivatives of the Standard Logistic Function

The `k`-th derivative of \\\sigma(z) = 1/(1 + e^{-z})\\, written as a
polynomial in \\p = \sigma(z)\\ itself.

## Usage

``` r
logistic_deriv(p, k)
```

## Arguments

- p:

  A numeric vector of logistic values, \\p = \sigma(z)\\.

- k:

  The derivative order, an integer from 1 to 4.

## Value

A numeric vector of the same length as `p`.

## Details

Three separate links need these same four polynomials, which is why they
live in one place rather than being written out three times:

- [`logit_link`](https://statmodels7.github.io/linkfunctions7/reference/logit_link.md)
  uses them directly, \\h^{(k)} = \sigma^{(k)}\\;

- [`bounded_link`](https://statmodels7.github.io/linkfunctions7/reference/bounded_link.md)
  with both endpoints scales them by the interval width, \\h^{(k)} = W
  \sigma^{(k)}\\;

- [`softplus_link`](https://statmodels7.github.io/linkfunctions7/reference/softplus_link.md)
  uses them shifted one order down, since the softplus is an
  antiderivative of the logistic: \\h^{(k+1)} = a^k \sigma^{(k)}\\.

The polynomials are \$\$\sigma' = p(1-p)\$\$ \$\$\sigma'' =
p(1-p)(1-2p)\$\$ \$\$\sigma''' = p(1-p)(1 - 6p + 6p^2)\$\$
\$\$\sigma'''' = p(1-p)(1 - 14p + 36p^2 - 24p^3)\$\$ and are evaluated
in Horner form, which is twice as fast at the fourth order and agrees
with the expanded form to within one unit in the last place.

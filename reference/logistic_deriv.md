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

  The derivative order, an integer from 1 to 5.

## Value

A numeric vector of the same length as `p`.

## Details

Three links use these polynomials:

- [`logit_link()`](https://statmodels7.github.io/linkfunctions7/reference/logit_link.md)
  uses them directly, \\h^{(k)} = \sigma^{(k)}\\;

- [`bounded_link()`](https://statmodels7.github.io/linkfunctions7/reference/bounded_link.md)
  with both endpoints scales them by the interval width, \\h^{(k)} = W
  \sigma^{(k)}\\;

- [`softplus_link()`](https://statmodels7.github.io/linkfunctions7/reference/softplus_link.md)
  uses them shifted one order down, since the softplus is an
  antiderivative of the logistic: \\h^{(k+1)} = a^k \sigma^{(k)}\\.

The three links evaluate these polynomials in compiled code, written in
\\p\\, \\q = 1 - p\\ and \\pq\\ with \\q\\ computed from \\z\\ and not
as \\1 - p\\, so that the derivatives keep their relative accuracy as
\\p\\ approaches 1: \\\sigma' = pq\\, \\\sigma'' = pq(q - p)\\,
\\\sigma''' = pq(1 - 6pq)\\, \\\sigma'''' = pq(q - p)(1 - 12pq)\\ and
\\\sigma^{(5)} = pq(1 - 30pq + 120p^2q^2)\\. The transcription of the
polynomials in \\p\\ alone is kept in `src/link_points.h` and reached
through `lk_logistic_poly_cpp()`. This function is the R statement of
the same five polynomials, and `test-logistic-twin.R` compares the two
at every order and checks that each of the three links reaches the
polynomial named in its description.

The comparison uses a tolerance instead of exact equality, because both
forms are Horner evaluations that contain multiply-adds, and a compiler
may fuse these into FMA instructions that skip an intermediate rounding.

The polynomials are \$\$\sigma' = p(1-p)\$\$ \$\$\sigma'' =
p(1-p)(1-2p)\$\$ \$\$\sigma''' = p(1-p)(1 - 6p + 6p^2)\$\$
\$\$\sigma'''' = p(1-p)(1 - 14p + 36p^2 - 24p^3)\$\$ \$\$\sigma^{(5)} =
p(1-p)(1 - 30p + 150p^2 - 240p^3 + 120p^4)\$\$ and are evaluated in
Horner form, which is faster than the expanded form and agrees with it
up to rounding error; near a root of a polynomial that error can be
large relative to the value. Each polynomial follows from the one before
by \\P\_{k+1} = (1-2p)P_k + p(1-p)P_k'\\, so a higher order can be
generated from the recurrence.

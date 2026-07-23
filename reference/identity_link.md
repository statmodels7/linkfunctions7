# The Identity Link Function

Creates an S7 object of class `link` implementing the Identity
transformation. This is the canonical link function for the mean
parameter of the Normal (Gaussian) distribution.

## Usage

``` r
identity_link()
```

## Value

An S7 object of class `IdentityLink` (inheriting from `link`) containing
the transformation functions and their exact analytical derivatives up
to the fourth order.

## Details

The Identity link is defined simply as \\\eta = \theta\\. Consequently,
the inverse link is also \\\theta = \eta\\.

All first derivatives are constant (equal to 1), and all higher-order
derivatives up to the fourth order are exactly zero.

The domain of \\\theta\\ is unbounded, meaning the valid domain is
`c(-Inf, Inf)`.

## See also

[`link`](https://statmodels7.github.io/linkfunctions7/reference/link.md)

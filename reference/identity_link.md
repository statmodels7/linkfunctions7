# The Identity Link Function

The identity link \\\eta = \theta\\, for a parameter that is already
unconstrained.

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

## Examples

``` r
lk <- identity_link()
lk
#> S7 Link Object: identity
#>   - Parameter domain (theta): (-Inf, Inf)

linkfun(lk, c(-1, 0, 1))
#> [1] -1  0  1
linkinv(lk, c(-1, 0, 1))
#> [1] -1  0  1

# the first derivative is 1 and every higher one is 0 ...
dlinkfun(lk, c(-1, 0, 1))
#> [1] 1 1 1
d2linkfun(lk, c(-1, 0, 1))
#> [1] 0 0 0

# ... but missingness is still propagated, not swallowed by the constant
dlinkfun(lk, c(1, NA))
#> [1]  1 NA
```

# The Identity Link Function

The identity link \\\eta = \theta\\, for a parameter that is already
unconstrained.

## Usage

``` r
identity_link()
```

## Value

An S7 object of class `IdentityLink`, inheriting from
[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
whose methods compute the link, its inverse and their derivatives to the
fifth order.

## Details

The identity link is \\\eta = \theta\\, and its inverse is \\\theta =
\eta\\.

The first derivative is 1 in both directions, and every higher
derivative is 0.

The domain of \\\theta\\ is the whole real line.

## See also

[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md)

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

# The General Bounded Link Function

The link for a parameter confined to \\(lwr, upr)\\: a scaled logit when
both endpoints are finite, a shifted log \\\eta = \log(\theta - lwr)\\
when only the lower is, its mirror image \\\eta = \log(upr - \theta)\\
when only the upper is, and the identity when neither is given.

## Usage

``` r
bounded_link(lwr = NULL, upr = NULL)
```

## Arguments

- lwr:

  Numeric or `NULL`. The lower bound of the interval.

- upr:

  Numeric or `NULL`. The upper bound of the interval.

## Value

An S7 object inheriting from
[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
whose methods compute the link, its inverse and their derivatives to the
fifth order. Its class depends on the endpoints given:
[`DoublyBoundedLink()`](https://statmodels7.github.io/linkfunctions7/reference/DoublyBoundedLink.md),
[`LowerBoundedLink()`](https://statmodels7.github.io/linkfunctions7/reference/LowerBoundedLink.md),
[`UpperBoundedLink()`](https://statmodels7.github.io/linkfunctions7/reference/UpperBoundedLink.md),
or
[`IdentityLink()`](https://statmodels7.github.io/linkfunctions7/reference/IdentityLink.md)
when neither endpoint is supplied.

## Details

The link depends on which endpoints are given:

- `lwr` and `upr`: \\\theta\\ is mapped to \\p = (\theta -
  \text{lwr})/(\text{upr} - \text{lwr})\\ in \\(0, 1)\\, and then
  through the logit;

- `lwr` only: \\\eta = \log(\theta - \text{lwr})\\, with inverse
  \\\theta = \exp(\eta) + \text{lwr}\\;

- `upr` only: \\\eta = \log(\text{upr} - \theta)\\, with inverse
  \\\theta = \text{upr} - \exp(\eta)\\;

- neither: the result is
  [`identity_link()`](https://statmodels7.github.io/linkfunctions7/reference/identity_link.md).

## See also

[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
[`logit_link()`](https://statmodels7.github.io/linkfunctions7/reference/logit_link.md),
[`log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.md)

## Examples

``` r
# a parameter known to lie in (0, 10)
lk <- bounded_link(lwr = 0, upr = 10)
lk
#> S7 Link Object: bounded(lwr=0, upr=10)
#>   - Parameter domain (theta): (0, 10)
#>   - Link parameters: lwr = 0, upr = 10
linkinv(lk, c(-2, 0, 2))     # always inside the interval
#> [1] 1.192029 5.000000 8.807971
linkfun(lk, c(1, 5, 9))
#> [1] -2.197225  0.000000  2.197225

# one-sided: a variance component bounded below by zero
bounded_link(lwr = 0)
#> S7 Link Object: lower_bounded(lwr=0)
#>   - Parameter domain (theta): (0, Inf)
#>   - Link parameters: lwr = 0

# no endpoints at all is the identity
bounded_link()
#> S7 Link Object: identity
#>   - Parameter domain (theta): (-Inf, Inf)

# derivatives, as for any other link
dlinkinv(lk, 0)
#> [1] 2.5
```

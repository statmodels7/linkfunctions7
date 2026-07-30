# The Power Link Function

Creates an S7 object of class `link` implementing the Power
transformation family. This function generates a specific link function
based on the provided power parameter `lambda`. It elegantly includes
the special case where `lambda = 0`, which dynamically returns the Log
link.

## Usage

``` r
power_link(lambda = 1)
```

## Arguments

- lambda:

  A numeric value defining the power of the transformation. Defaults to
  1.

## Value

An S7 object of class `PowerLink` (inheriting from `link`), or an object
of class `LogLink` if `lambda = 0`.

## Details

The Power link is defined mathematically as \\\eta = \theta^\lambda\\.
Consequently, the inverse link is derived as \\\theta =
\eta^{1/\lambda}\\.

**Special Case (Box-Cox continuity):** If `lambda = 0`, the function
mathematically approaches \\\log(\theta)\\. In this scenario, the
function automatically instantiates and returns a
[`log_link`](https://statmodels7.github.io/linkfunctions7/reference/log_link.md)
object, modifying its internal state to reflect the `lambda = 0`
parameter.

Common special cases include:

- `lambda = 1`: Identity link.

- `lambda = 0.5`: Square-root link.

- `lambda = -1`: Inverse link.

- `lambda = 0`: Log link.

The mathematical domain of \\\theta\\ is `c(0, Inf)`. Depending on the
value of `lambda`, extreme care must be taken during numerical
optimization to guarantee that \\\eta\\ remains strictly positive to
avoid `NaN`s from fractional exponents.

## See also

[`link`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
[`log_link`](https://statmodels7.github.io/linkfunctions7/reference/log_link.md),
[`identity_link`](https://statmodels7.github.io/linkfunctions7/reference/identity_link.md)

## Examples

``` r
lk <- power_link(2)
lk
#> S7 Link Object: power(lambda=2)
#>   - Parameter domain (theta): (0, Inf)
#>   - Link parameters: lambda = 2

theta <- c(1, 2, 3)
eta <- linkfun(lk, theta)
eta
#> [1] 1 4 9
linkinv(lk, eta)
#> [1] 1 2 3

# special cases of the family
linkfun(power_link(1),    5)   # identity
#> [1] 5
linkfun(power_link(0.5),  4)   # square root
#> [1] 2
linkfun(power_link(-1),   4)   # inverse
#> [1] 0.25

# lambda = 0 is the log link by continuity, and is returned as one
power_link(0)
#> S7 Link Object: power(lambda=0)
#>   - Parameter domain (theta): (0, Inf)
#>   - Link parameters: lambda = 0
linkfun(power_link(0), exp(1))
#> [1] 1
```

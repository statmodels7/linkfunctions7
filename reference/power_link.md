# The Power Link Function

The power link \\\eta = \theta^\lambda\\ on \\(0, \infty)\\; at
\\\lambda = 0\\ it returns the log link, the limit of the Box-Cox
transformation.

## Usage

``` r
power_link(lambda = 1)
```

## Arguments

- lambda:

  A numeric value defining the power of the transformation. Defaults to
  1.

## Value

An S7 object of class `PowerLink`, inheriting from
[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
or of class `LogLink` when `lambda = 0`.

## Details

The power link is \\\eta = \theta^\lambda\\, with inverse \\\theta =
\eta^{1/\lambda}\\.

As \\\lambda \to 0\\ the Box-Cox transformation \\(\theta^\lambda -
1)/\lambda\\ tends to \\\log(\theta)\\, so for `lambda = 0` the function
returns a
[`log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.md)
object whose `link_params` record `lambda = 0`.

Common special cases are

- `lambda = 1`: identity link;

- `lambda = 0.5`: square-root link;

- `lambda = -1`: inverse link;

- `lambda = 0`: log link.

The domain of \\\theta\\ is \\(0, \infty)\\. For every non-zero `lambda`
the image of the link is also \\(0, \infty)\\, so the linear predictor
must stay positive during optimization;
[`eta_bounds()`](https://statmodels7.github.io/linkfunctions7/reference/eta_bounds.md)
returns this range.

## See also

[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
[`log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.md),
[`identity_link()`](https://statmodels7.github.io/linkfunctions7/reference/identity_link.md)

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

# lambda = 0 returns the log link, the limit of the Box-Cox transformation
power_link(0)
#> S7 Link Object: power(lambda=0)
#>   - Parameter domain (theta): (0, Inf)
#>   - Link parameters: lambda = 0
linkfun(power_link(0), exp(1))
#> [1] 1
```

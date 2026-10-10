# The Probit Link Function

The probit link \\\eta = \Phi^{-1}(\theta)\\ on \\(0, 1)\\, with
\\\Phi\\ the standard normal distribution function.

## Usage

``` r
probit_link()
```

## Value

An S7 object of class `ProbitLink`, inheriting from
[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
whose methods compute the link, its inverse and their derivatives to the
fifth order.

## Details

The probit link is \\\eta = \Phi^{-1}(\theta)\\, where \\\Phi^{-1}\\ is
the quantile function of the standard normal distribution
([`qnorm()`](https://rdrr.io/r/stats/Normal.html)). The inverse link is
the standard normal distribution function, \\\theta = \Phi(\eta)\\
([`pnorm()`](https://rdrr.io/r/stats/Normal.html)).

Like the logit, the probit link is symmetric about \\\theta = 1/2\\,
where \\\eta = 0\\. The tails of the normal distribution are lighter
than those of the logistic, so \\\theta\\ approaches 0 and 1 faster than
under the logit link.

The domain of \\\theta\\ is \\(0, 1)\\.

## See also

[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
[`logit_link()`](https://statmodels7.github.io/linkfunctions7/reference/logit_link.md),
[`cauchit_link()`](https://statmodels7.github.io/linkfunctions7/reference/cauchit_link.md)

## Examples

``` r
lk <- probit_link()
lk
#> S7 Link Object: probit
#>   - Parameter domain (theta): (0, 1)

p <- c(0.1, 0.5, 0.9)
eta <- linkfun(lk, p)      # standard normal quantiles
eta
#> [1] -1.281552  0.000000  1.281552
linkinv(lk, eta)
#> [1] 0.1 0.5 0.9

# the first inverse derivative is the standard normal density
dlinkinv(lk, 0)
#> [1] 0.3989423
dnorm(0)
#> [1] 0.3989423

# probit tails approach 0 and 1 faster than logit ones
linkinv(probit_link(), 3)
#> [1] 0.9986501
linkinv(logit_link(), 3)
#> [1] 0.9525741
```

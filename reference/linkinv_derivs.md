# Derivatives of an Inverse Link up to an Order

Computes the derivatives of the inverse link \\g^{-1}(\eta)\\ of orders
one to `order` at the same linear predictors, as a list.

## Usage

``` r
linkinv_derivs(x, eta, order)
```

## Arguments

- x:

  An object of class `link`.

- eta:

  A numeric vector of linear predictors, as for
  [`dlinkinv()`](https://statmodels7.github.io/linkfunctions7/reference/dlinkinv.md).

- order:

  The highest derivative order, a single whole number from 1 to 6.

## Value

A list of length `order` whose element `k` is the `k`-th derivative of
\\g^{-1}\\ at `eta`, a numeric vector of the same length as `eta`.

## Details

The method of the base class calls
[`dlinkinv()`](https://statmodels7.github.io/linkfunctions7/reference/dlinkinv.md),
[`d2linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/d2linkinv.md)
and the higher generics in turn, so the numerical fallbacks of a link
apply here as they do there. A class whose derivatives share their work
may register its own method: the log link computes \\e^\eta\\ once and
returns that vector at every order, which gives the same values as the
six generics at a sixth of the cost. A class that inherits from the log
link and overrides one of its derivatives is evaluated through the
generics.

## See also

[`dlinkinv()`](https://statmodels7.github.io/linkfunctions7/reference/dlinkinv.md)
and its siblings for one order at a time,
[`linkinvderiv()`](https://statmodels7.github.io/linkfunctions7/reference/linkinvderiv.md)
for one order chosen at run time.

## Examples

``` r
# the four derivatives of exp at 0 are all equal to 1
linkinv_derivs(log_link(), 0, order = 4)
#> [[1]]
#> [1] 1
#> 
#> [[2]]
#> [1] 1
#> 
#> [[3]]
#> [1] 1
#> 
#> [[4]]
#> [1] 1
#> 

# the logistic's first two derivatives at eta = 0
linkinv_derivs(logit_link(), 0, order = 2)
#> [[1]]
#> [1] 0.25
#> 
#> [[2]]
#> [1] 0
#> 
```

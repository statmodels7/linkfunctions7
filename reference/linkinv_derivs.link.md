# Inverse Link Derivatives from the Order-Specific Generics

Evaluates
[`dlinkinv()`](https://statmodels7.github.io/linkfunctions7/reference/dlinkinv.md)
and the higher generics up to `order` and returns them as a list. It is
the method of
[`linkinv_derivs()`](https://statmodels7.github.io/linkfunctions7/reference/linkinv_derivs.md)
for the base class.

## Usage

``` r
linkinv_derivs.link(x, eta, order)
```

## Arguments

- x:

  An object of class `link`.

- eta:

  A numeric vector.

- order:

  A whole number from 1 to 5.

## Value

A list of `order` numeric vectors of the same length as `eta`.

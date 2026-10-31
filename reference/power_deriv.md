# Derivative of a Power

Computes the derivative of order `order` of \\v^p\\, which is
\\c\\v^{p - \mathrm{order}}\\ with \\c = p(p-1)\cdots(p -
\mathrm{order} + 1)\\. It is the body of the derivative methods of
[`power_link()`](https://statmodels7.github.io/linkfunctions7/reference/power_link.md),
in both directions.

## Usage

``` r
power_deriv(v, p, order)
```

## Arguments

- v:

  A numeric vector: \\\theta\\ in the forward direction, \\\eta\\ in the
  inverse one.

- p:

  The exponent: \\\lambda\\ forward, \\1/\lambda\\ inverse.

- order:

  The derivative order, from 1 to 6.

## Value

A numeric vector of the same length as `v`, missing wherever `v` is.

## Details

The coefficient \\c\\ is multiplied from left to right. When it is
exactly zero, which happens for a whole number \\p\\ from 0 to
`order - 1`, the derivative is identically zero and zeros are returned
directly, so that the product \\0 \cdot v^{p - \mathrm{order}}\\ does
not give `NaN` where the power overflows. A product that overflows or
underflows in its second factor alone is recomputed with the power split
in two halves, \\(c\\v^h)\\v^h\\ with \\h = (p - \mathrm{order})/2\\,
which is finite whenever the value is.

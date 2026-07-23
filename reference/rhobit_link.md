# The Rhobit (Fisher's z) Link Function

Creates an S7 object of class `link` implementing the Rhobit
transformation, also known as Fisher's z-transformation. This link
function rigorously maps the open interval `c(-1, 1)` to the real line
`c(-Inf, Inf)`. It is primarily used for modeling correlation
coefficients or other bounded parameters that are symmetrically
constrained.

## Usage

``` r
rhobit_link()
```

## Value

An S7 object of class `RhobitLink` (inheriting from `link`) containing
the transformation functions and their exact analytical derivatives up
to the fourth order.

## Details

The Rhobit link is defined mathematically using the inverse hyperbolic
tangent function: \\\eta = \text{arctanh}(\theta) = \frac{1}{2}
\log\left(\frac{1 + \theta}{1 - \theta}\right)\\.

The inverse link is the hyperbolic tangent function: \\\theta =
\tanh(\eta) = \frac{\exp(2\eta) - 1}{\exp(2\eta) + 1}\\.

The valid mathematical domain of \\\theta\\ is exactly `c(-1, 1)`.

## See also

[`link`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
[`logit_link`](https://statmodels7.github.io/linkfunctions7/reference/logit_link.md)

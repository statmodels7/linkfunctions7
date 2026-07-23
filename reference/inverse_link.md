# The Inverse (Reciprocal) Link Function

Creates an S7 object of class `link` implementing the reciprocal
transformation.

## Usage

``` r
inverse_link()
```

## Value

An S7 object of class `InverseLink` (inheriting from `link`) containing
the transformation functions and their exact analytical derivatives up
to the fourth order.

## Details

The Inverse link is defined as \\\eta = 1/\theta\\. The inverse link
function is therefore perfectly symmetric: \\\theta = 1/\eta\\.

This link is typically used for modeling positive continuous data where
the mean is inversely proportional to the linear predictor (e.g., in
Gamma regression).

The domain of \\\theta\\ is conventionally `c(0, Inf)`. Care must be
taken to ensure the linear predictor \\\eta\\ remains strictly positive
(or strictly negative) during optimization to avoid division by zero or
mapping to invalid negative parameter values.

## See also

[`link`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
[`identity_link`](https://statmodels7.github.io/linkfunctions7/reference/identity_link.md)

# The Inverse Square Link Function

Creates an S7 object of class `link` implementing the Inverse Square
transformation. This is the canonical link function for the Inverse
Gaussian distribution.

## Usage

``` r
inverse_sq_link()
```

## Value

An S7 object of class `InverseSqLink` (inheriting from `link`)
containing the transformation functions and their exact analytical
derivatives up to the fourth order.

## Details

The Inverse Square link is defined mathematically as \\\eta = 1 /
\theta^2\\. Consequently, the inverse link function is derived as
\\\theta = 1 / \sqrt{\eta}\\.

This specific link function is predominantly utilized in Generalized
Linear Models (GLMs) assuming an Inverse Gaussian response distribution.
In such frameworks, the variance is proportional to the cube of the mean
(\\\text{Var}(Y) \propto \mu^3\\).

**Domain and Optimization Constraints:** Both the parameter \\\theta\\
and the linear predictor \\\eta\\ must be strictly positive. The valid
mathematical domain for \\\theta\\ is `c(0, Inf)`. During optimization
routines (e.g., Fisher Scoring or Newton-Raphson), extreme care must be
taken to ensure the linear predictor \\\eta \> 0\\. Evaluating the
inverse link or its derivatives at non-positive values of \\\eta\\ will
inevitably result in `NaN`s due to fractional powers and square root
operations.

## See also

[`link`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
[`inverse_link`](https://statmodels7.github.io/linkfunctions7/reference/inverse_link.md)

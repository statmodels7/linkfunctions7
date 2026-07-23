# The Softplus Link Function

Creates an S7 object of class `link` implementing the Softplus
transformation. The Softplus function is a smooth approximation of the
rectifier (ReLU) and ensures the parameter \\\theta\\ remains positive.
Unlike the Log link, which implies an exponential relationship
throughout, the Softplus link asymptotically approaches linearity for
large positive values of the linear predictor.

## Usage

``` r
softplus_link(a = 1)
```

## Arguments

- a:

  A numeric value specifying the scaling parameter
  (smoothness/steepness). Must be strictly positive. Defaults to 1.

## Value

An S7 object of class `SoftplusLink` (inheriting from `link`) containing
the transformation functions, their exact analytical derivatives up to
the fourth order, and the parameter `a`.

## Details

The Softplus link describes the relationship where the response
parameter \\\theta\\ is the Softplus of the linear predictor \\\eta\\.

Mathematically:

- Inverse Link (Softplus): \\\theta = \frac{1}{a} \log(1 + \exp(a
  \eta))\\

- Link Function: \\\eta = \frac{1}{a} \log(\exp(a \theta) - 1)\\

**Behavior:** For large negative \\\eta\\, \\\theta \approx 0\\. For
large positive \\\eta\\, \\\theta \approx \eta\\ (linear behavior),
whereas a Log link would imply \\\theta = \exp(\eta)\\ (exponential
behavior).

**Numerical Stability:** The inverse link implementation intelligently
utilizes conditional algebraic logic to ensure robust numerical
stability for large positive values of \\\eta\\, entirely avoiding
precision overflow.

The mathematical domain of \\\theta\\ is `c(0, Inf)`.

## See also

[`link`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
[`log_link`](https://statmodels7.github.io/linkfunctions7/reference/log_link.md),
[`identity_link`](https://statmodels7.github.io/linkfunctions7/reference/identity_link.md)

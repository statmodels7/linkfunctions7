# The Square Root Link Function

Creates an S7 object of class `link` implementing the Square Root
transformation. This is a specific case of the Power link family (with
\\\lambda = 0.5\\) and is prominently used for modeling count data
(e.g., Poisson regression) as a variance-stabilizing transformation.

## Usage

``` r
sqrt_link()
```

## Value

An S7 object of class `SqrtLink` (inheriting from `link`) containing the
transformation functions and their exact analytical derivatives up to
the fourth order.

## Details

The Square Root link is mathematically defined as \\\eta =
\sqrt{\theta}\\. Consequently, the inverse link is derived as \\\theta =
\eta^2\\.

Unlike the Log link, this transformation allows \\\theta\\ to reach 0
exactly. While the inverse function (\\\eta^2\\) is mathematically valid
for negative values of \\\eta\\, in the specific context of this link
function, the linear predictor \\\eta\\ is typically constrained to be
non-negative. This restriction preserves a strictly one-to-one mapping
with \\\theta\\.

The strict mathematical domain for \\\theta\\ is `c(0, Inf)`.

## See also

[`link`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
[`power_link`](https://statmodels7.github.io/linkfunctions7/reference/power_link.md),
[`log_link`](https://statmodels7.github.io/linkfunctions7/reference/log_link.md)

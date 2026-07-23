# The Probit Link Function

Creates an S7 object of class `link` implementing the Probit
transformation. This link function rigorously relies on the cumulative
distribution function (CDF) of the standard Normal distribution. It is
widely used in Generalized Linear Models (GLMs) for binary data, often
justified by a latent normal variable interpretation.

## Usage

``` r
probit_link()
```

## Value

An S7 object of class `ProbitLink` (inheriting from `link`) containing
the transformation functions and their exact analytical derivatives up
to the fourth order.

## Details

The Probit link is mathematically defined as \\\eta =
\Phi^{-1}(\theta)\\, where \\\Phi^{-1}\\ is the quantile function of the
standard normal distribution (`qnorm`). The inverse link is \\\theta =
\Phi(\eta)\\, the standard normal CDF (`pnorm`).

Similarly to the `logit` link, the Probit is symmetric around \\\theta =
0.5\\ (where \\\eta = 0\\). However, the tails of the Normal
distribution approach 0 and 1 faster than the Logistic distribution.

The strictly mathematical domain of \\\theta\\ is `c(0, 1)`.

## See also

[`link`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
[`logit_link`](https://statmodels7.github.io/linkfunctions7/reference/logit_link.md),
[`cauchit_link`](https://statmodels7.github.io/linkfunctions7/reference/cauchit_link.md)

# The Logit Link Function

Creates an S7 object of class `link` implementing the Logit (log-odds)
transformation. This is the canonical link function for the success
probability parameter of the Bernoulli and Binomial distributions and
serves as the foundation of Logistic Regression.

## Usage

``` r
logit_link()
```

## Value

An S7 object of class `LogitLink` (inheriting from `link`) containing
the transformation functions and their exact analytical derivatives up
to the fourth order.

## Details

The Logit link is defined mathematically as \\\eta =
\log(\frac{\theta}{1 - \theta})\\. The inverse link is the standard
logistic function (sigmoid): \\\theta = \frac{1}{1 + \exp(-\eta)}\\.

The link is perfectly symmetric around \\\theta = 0.5\\ (which
corresponds to \\\eta = 0\\). It elegantly interprets the linear
predictor as the log-odds of the event probability.

The mathematical domain of \\\theta\\ is `c(0, 1)`.

**Implementation Details:** This function internally delegates to R's
highly optimized native functions
[`stats::qlogis`](https://rdrr.io/r/stats/Logistic.html) and
[`stats::plogis`](https://rdrr.io/r/stats/Logistic.html). This ensures
maximum numerical stability and precision, especially when evaluating
probabilities exceedingly close to the boundaries of 0 and 1.

## See also

[`link`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
[`probit_link`](https://statmodels7.github.io/linkfunctions7/reference/probit_link.md),
[`cloglog_link`](https://statmodels7.github.io/linkfunctions7/reference/cloglog_link.md)

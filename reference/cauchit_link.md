# The Cauchit Link Function

Creates an S7 object of class `link` implementing the Cauchit
transformation. This link function rigorously maps the open interval
`c(0, 1)` to the real line by utilizing the quantile function of the
standard Cauchy distribution.

## Usage

``` r
cauchit_link()
```

## Value

An S7 object of class `CauchitLink` (inheriting from `link`) containing
the transformation functions and their exact analytical derivatives up
to the fourth order.

## Details

The Cauchit link is defined mathematically as \\\eta = \tan(\pi(\theta -
0.5))\\, which corresponds perfectly to `qcauchy(theta)`. The inverse
link is the standard Cauchy CDF \\\theta = \frac{1}{\pi} \arctan(\eta) +
0.5\\, computed via `pcauchy(eta)`.

**Heavy Tails:** Unlike the Logit or Probit links, the Cauchit link has
exceedingly heavier tails. This makes it particularly robust and useful
for modeling binary data where the probability approaches 0 or 1 very
slowly, or when the dataset contains severe outliers that might
disproportionately influence the fit of light-tailed link functions.

The strictly valid mathematical domain for \\\theta\\ is `c(0, 1)`.

## See also

[`link`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
[`logit_link`](https://statmodels7.github.io/linkfunctions7/reference/logit_link.md),
[`probit_link`](https://statmodels7.github.io/linkfunctions7/reference/probit_link.md)

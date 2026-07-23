# The Complementary Log-Log (ClogLog) Link Function

Creates an S7 object of class `link` implementing the Complementary
Log-Log transformation. This link is highly asymmetric and is
predominantly used for modeling binary data where the probability of the
event approaches 1 very slowly but approaches 0 rather sharply.

## Usage

``` r
cloglog_link()
```

## Value

An S7 object of class `ClogLogLink` (inheriting from `link`) and their
exact analytical derivatives up to the fourth order.

## Details

The ClogLog link is defined mathematically as \\\eta = \log(-\log(1 -
\theta))\\. Consequently, the inverse link is derived as \\\theta = 1 -
\exp(-\exp(\eta))\\.

Unlike the symmetric Logit and Probit links, the ClogLog link lacks
symmetry. It is fundamentally related to the Extreme Value (Gumbel)
distribution and is frequently utilized in discrete-time survival
analysis (proportional hazards models) as well as for modeling rare
events.

The strictly valid mathematical domain for \\\theta\\ is `c(0, 1)`.

## See also

[`link`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
[`logit_link`](https://statmodels7.github.io/linkfunctions7/reference/logit_link.md),
[`loglog_link`](https://statmodels7.github.io/linkfunctions7/reference/loglog_link.md)

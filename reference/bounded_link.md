# The General Bounded Link Function

Creates an S7 object of class `link` that maps a constrained interval to
the whole real line. By specifying `lwr` and `upr`, this function
dynamically constructs a doubly bounded (scaled logit), lower bounded
(shifted log), upper bounded, or unbounded (identity) link function.

## Usage

``` r
bounded_link(lwr = NULL, upr = NULL)
```

## Arguments

- lwr:

  Numeric or `NULL`. The lower bound of the interval.

- upr:

  Numeric or `NULL`. The upper bound of the interval.

## Value

An S7 object of class `link` containing the transformation functions and
their exact analytical derivatives up to the fourth order.

## Details

**Doubly Bounded (`lwr` and `upr` provided):** Transforms \\\theta\\ by
normalizing it to `c(0, 1)` via \\p = \frac{\theta -
\text{lwr}}{\text{upr} - \text{lwr}}\\, and then applying the logit
function.

**Lower Bounded (`lwr` provided, `upr = NULL`):** Defined as \\\eta =
\log(\theta - \text{lwr})\\, with inverse \\\theta = \exp(\eta) +
\text{lwr}\\.

**Upper Bounded (`lwr = NULL`, `upr` provided):** Defined as \\\eta =
\log(\text{upr} - \theta)\\, with inverse \\\theta = \text{upr} -
\exp(\eta)\\.

**Unbounded (`lwr = NULL`, `upr = NULL`):** Returns the standard
[`identity_link`](https://statmodels7.github.io/linkfunctions7/reference/identity_link.md).

## See also

[`link`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
[`logit_link`](https://statmodels7.github.io/linkfunctions7/reference/logit_link.md),
[`log_link`](https://statmodels7.github.io/linkfunctions7/reference/log_link.md)

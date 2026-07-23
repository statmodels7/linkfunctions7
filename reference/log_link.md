# The Logarithmic Link Function

Creates an S7 object of class `link` implementing the natural logarithm
transformation. This is the canonical link function for the mean of the
Poisson distribution and is widely used for modeling count data or
non-negative continuous data with multiplicative effects.

## Usage

``` r
log_link()
```

## Value

An S7 object of class `LogLink` (inheriting from `link`) containing the
transformation functions and their exact analytical derivatives up to
the fourth order.

## Details

The Log link is defined mathematically as \\\eta = \log(\theta)\\. The
inverse link is the exponential function \\\theta = \exp(\eta)\\.

A remarkable mathematical property of this link is that the inverse
function is its own derivative. Therefore, the parameter \\\theta\\ and
all its derivatives with respect to \\\eta\\ are equal to
\\\exp(\eta)\\.

The valid mathematical domain of \\\theta\\ is `c(0, Inf)`.

**Numerical Stability:** During the evaluation of the inverse link and
its derivatives, the result is bounded from below by
`.Machine$double.eps`. This prevents numerical underflow to exactly zero
when \\\eta\\ is a large negative number, which would otherwise produce
`Inf` when subsequently calculating \\1/\theta\\.

## See also

[`link`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
[`inverse_link`](https://statmodels7.github.io/linkfunctions7/reference/inverse_link.md)

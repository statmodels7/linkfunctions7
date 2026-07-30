# S7 Class for Statistical Link Functions

A strictly typed S7 object that encapsulates the metadata of a
statistical link function. The mathematical transformations, including
the forward and inverse functions and their exact analytical derivatives
up to the fourth order, are implemented and registered as S7 generic
methods.

## Usage

``` r
link(link_name = character(0), link_bounds = integer(0), link_params = NULL)
```

## Arguments

- link_name:

  A character string identifying the link (e.g., "logit").

- link_bounds:

  A numeric vector of length 2 `c(lower, upper)` defining the valid
  domain for \\\theta\\.

- link_params:

  A list or vector of additional parameters required to define the link,
  or `NULL`.

## Value

An S7 object of class `link`. In practice this class is not instantiated
directly: each link is a subclass created by one of the constructors
([`logit_link`](https://statmodels7.github.io/linkfunctions7/reference/logit_link.md),
[`power_link`](https://statmodels7.github.io/linkfunctions7/reference/power_link.md),
...), and `link` is what they all inherit from and what methods dispatch
on.

## Details

Objects of class `link` are instantiated using the S7 object system.

The object assumes the following mathematical notation:

- \\\theta\\: The response parameter (e.g., probability, mean,
  dispersion).

- \\\eta\\: The linear predictor (unconstrained scale).

The relationship is defined as \\\eta = g(\theta)\\ (link function) and
\\\theta = g^{-1}(\eta)\\ (inverse link function).

## Examples

``` r
# every constructor returns an object inheriting from `link`
lk <- logit_link()
lk
#> S7 Link Object: logit
#>   - Parameter domain (theta): (0, 1)
S7::S7_inherits(lk, link)
#> [1] TRUE

lk@link_name
#> [1] "logit"
lk@link_bounds
#> [1] 0 1
```

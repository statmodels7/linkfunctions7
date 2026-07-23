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

## Details

Objects of class `link` are instantiated using the S7 object system.

The object assumes the following mathematical notation:

- \\\theta\\: The response parameter (e.g., probability, mean,
  dispersion).

- \\\eta\\: The linear predictor (unconstrained scale).

The relationship is defined as \\\eta = g(\theta)\\ (link function) and
\\\theta = g^{-1}(\eta)\\ (inverse link function).

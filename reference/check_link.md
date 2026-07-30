# Validate and Check a Link Object

A diagnostic S7 method to mathematically validate a `link` object. It
sequentially verifies the algebraic invertibility and the correctness of
the analytical derivatives using numerical gradients in a chained
sequence.

## Usage

``` r
check_link(x, tolerance = 1e-05, ...)

check_link.link(x, tolerance = 1e-05, ...)
```

## Arguments

- x:

  An object of class `link`.

- tolerance:

  Numeric tolerance for floating-point comparisons.

- ...:

  Additional arguments passed to methods.

## Value

Invisibly, a named list of check results; see `check_link.link` for its
shape. Called mainly for the summary printed to the console.

Invisibly, a named list of the check results: the four scalar logicals
`invertibility_theta`, `invertibility_eta`, `monotonicity` and
`inverse_theorem`, plus `link_derivatives` and
`inverse_link_derivatives`, each a logical vector of length four named
`order_1` to `order_4`. In those two, `TRUE` and `FALSE` mean what they
say and `NA` means **not checked**: the order is supplied by a numerical
fallback, so the value and the reference would be the same arithmetic
and would agree whatever the link did. The number of orders actually
implemented is carried on the result as the attribute
`"analytic_orders"`; see
[`link_fallback_orders`](https://statmodels7.github.io/linkfunctions7/reference/link_fallback_orders.md).
A derivative that raises an error still counts as `FALSE`. Called mainly
for the summary printed to the console.

## Details

The function assumes the existence of S7 generics `linkfun`, `linkinv`,
`linkderiv`, and `linkinvderiv`. The method performs the following six
diagnostic checks:

1.  **Invertibility (\\\theta\\ space):** Verifies \\g^{-1}(g(\theta)) =
    \theta\\. Ensures that mapping from the parameter space to the
    linear predictor and back is lossless.

2.  **Invertibility (\\\eta\\ space):** Verifies \\g(g^{-1}(\eta)) =
    \eta\\. Ensures that mapping from the linear predictor to the
    parameter space and back is lossless. Note that this test may fail
    intentionally and correctly for links that map to a restricted
    \\\eta\\ domain (e.g., the square root link).

3.  **Strict Monotonicity:** Checks if the first derivative
    \\g'(\theta)\\ is strictly positive or strictly negative across the
    domain, guaranteeing a one-to-one mapping.

4.  **Inverse Function Theorem:** Verifies the mathematical identity
    \\g'(\theta) \cdot (g^{-1})'(\eta) = 1\\, confirming the theoretical
    relationship between the link derivative and the inverse link
    derivative.

5.  **Link Derivatives:** Validates the exact analytical forward
    derivatives of \\g(\theta)\\ up to the 4th order by comparing them
    against numerical gradients.

6.  **Inverse Link Derivatives:** Validates the exact analytical inverse
    derivatives of \\g^{-1}(\eta)\\ up to the 4th order by comparing
    them against numerical gradients.

Both forward and inverse derivative testing avoids compounding numerical
errors by applying first-order numerical differentiation iteratively to
the exact lower-order analytical derivatives.

## Examples

``` r
check_link(sqrt_link())
#> Checking S7 Link Object: sqrt 
#>   [1] Invertibility (Theta space): [PASSED] 
#>   [2] Invertibility (Eta space):   [PASSED] 
#>   [3] Strict Monotonicity:         [PASSED] 
#>   [4] Inverse Function Theorem:    [PASSED] 
#>   [5] Link Derivatives:            [PASSED] 
#>   [6] Inverse Link Derivatives:    [PASSED] 
# every link the package ships passes all six checks
check_link(logit_link())
#> Checking S7 Link Object: logit 
#>   [1] Invertibility (Theta space): [PASSED] 
#>   [2] Invertibility (Eta space):   [PASSED] 
#>   [3] Strict Monotonicity:         [PASSED] 
#>   [4] Inverse Function Theorem:    [PASSED] 
#>   [5] Link Derivatives:            [PASSED] 
#>   [6] Inverse Link Derivatives:    [PASSED] 

res <- check_link(power_link(2))
#> Checking S7 Link Object: power(lambda=2) 
#>   [1] Invertibility (Theta space): [PASSED] 
#>   [2] Invertibility (Eta space):   [PASSED] 
#>   [3] Strict Monotonicity:         [PASSED] 
#>   [4] Inverse Function Theorem:    [PASSED] 
#>   [5] Link Derivatives:            [PASSED] 
#>   [6] Inverse Link Derivatives:    [PASSED] 
res$link_derivatives
#> order_1 order_2 order_3 order_4 
#>    TRUE    TRUE    TRUE    TRUE 
res$inverse_theorem
#> [1] TRUE

# the checks are what a user-defined link should be held to as well
all(unlist(check_link(bounded_link(0, 10))))
#> Checking S7 Link Object: bounded(lwr=0, upr=10) 
#>   [1] Invertibility (Theta space): [PASSED] 
#>   [2] Invertibility (Eta space):   [PASSED] 
#>   [3] Strict Monotonicity:         [PASSED] 
#>   [4] Inverse Function Theorem:    [PASSED] 
#>   [5] Link Derivatives:            [PASSED] 
#>   [6] Inverse Link Derivatives:    [PASSED] 
#> [1] TRUE
```

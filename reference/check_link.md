# Validate and Check a Link Object

Validates a `link` object numerically: invertibility in both directions
on a grid, strict monotonicity, the inverse function theorem
\\h'(\eta)\\g'(\theta) = 1\\, and every analytic derivative against one
numerical differentiation of the analytic order below it.

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

Invisibly, a named list of check results; see `check_link.link()` for
its shape. Called mainly for the summary printed to the console.

Invisibly, a named list of the check results: the four scalar logicals
`invertibility_theta`, `invertibility_eta`, `monotonicity` and
`inverse_theorem`, plus `link_derivatives` and
`inverse_link_derivatives`, each a logical vector of length five named
`order_1` to `order_5`. In those two, `NA` means that the order was not
checked: it is supplied by a numerical fallback, so the value and the
reference would come from the same computation and would agree even for
a wrong link. The number of orders actually implemented is carried on
the result as the attribute `"analytic_orders"`; see
[`link_fallback_orders()`](https://statmodels7.github.io/linkfunctions7/reference/link_fallback_orders.md).
A derivative that raises an error still counts as `FALSE`. Called mainly
for the summary printed to the console.

## Details

The method performs six checks:

1.  **Invertibility (\\\theta\\ space):** verifies \\g^{-1}(g(\theta)) =
    \theta\\ on a grid of parameter values.

2.  **Invertibility (\\\eta\\ space):** verifies \\g(g^{-1}(\eta)) =
    \eta\\ on a grid spanning the linear predictors that the link
    produces from the parameter grid.

3.  **Strict monotonicity:** checks that \\g'(\theta)\\ has the same
    sign at every point of the grid, so that the map is one-to-one.

4.  **Inverse function theorem:** verifies \\g'(\theta) \cdot
    (g^{-1})'(\eta) = 1\\.

5.  **Link derivatives:** compares each analytic derivative of
    \\g(\theta)\\, up to the fifth order, with a numerical one.

6.  **Inverse link derivatives:** compares each analytic derivative of
    \\g^{-1}(\eta)\\, up to the fifth order, with a numerical one.

Each analytic derivative of order \\k\\ is compared with one numerical
differentiation of the analytic derivative of order \\k - 1\\, so the
errors of successive numerical differentiations do not compound.

## See also

[`linkfun()`](https://statmodels7.github.io/linkfunctions7/reference/linkfun.md),
[`linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/linkinv.md),
[`linkderiv()`](https://statmodels7.github.io/linkfunctions7/reference/linkderiv.md),
[`linkinvderiv()`](https://statmodels7.github.io/linkfunctions7/reference/linkinvderiv.md),
[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md)

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
#> order_1 order_2 order_3 order_4 order_5 
#>    TRUE    TRUE    TRUE    TRUE    TRUE 
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

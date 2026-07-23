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

A logical list returning the success status of all available checks.

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

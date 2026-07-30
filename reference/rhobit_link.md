# The Rhobit (Fisher's z) Link Function

Creates an S7 object of class `link` implementing the Rhobit
transformation, also known as Fisher's z-transformation. This link
function rigorously maps the open interval `c(-1, 1)` to the real line
`c(-Inf, Inf)`. It is primarily used for modeling correlation
coefficients or other bounded parameters that are symmetrically
constrained.

## Usage

``` r
rhobit_link()
```

## Value

An S7 object of class `RhobitLink` (inheriting from `link`) containing
the transformation functions and their exact analytical derivatives up
to the fourth order.

## Details

The Rhobit link is defined mathematically using the inverse hyperbolic
tangent function: \\\eta = \text{arctanh}(\theta) = \frac{1}{2}
\log\left(\frac{1 + \theta}{1 - \theta}\right)\\.

The inverse link is the hyperbolic tangent function: \\\theta =
\tanh(\eta) = \frac{\exp(2\eta) - 1}{\exp(2\eta) + 1}\\.

The valid mathematical domain of \\\theta\\ is exactly `c(-1, 1)`.

## See also

[`link`](https://statmodels7.github.io/linkfunctions7/reference/link.md),
[`logit_link`](https://statmodels7.github.io/linkfunctions7/reference/logit_link.md)

## Examples

``` r
lk <- rhobit_link()
lk
#> S7 Link Object: rhobit
#>   - Parameter domain (theta): (-1, 1)

# built for parameters constrained to (-1, 1), such as a correlation
rho <- c(-0.9, 0, 0.9)
eta <- linkfun(lk, rho)    # Fisher's z
eta
#> [1] -1.472219  0.000000  1.472219
linkinv(lk, eta)
#> [1] -0.9  0.0  0.9

# the inverse is tanh, so its first derivative is the squared sech
dlinkinv(lk, 0)
#> [1] 1

check_link(lk)
#> Checking S7 Link Object: rhobit 
#>   [1] Invertibility (Theta space): [PASSED] 
#>   [2] Invertibility (Eta space):   [PASSED] 
#>   [3] Strict Monotonicity:         [PASSED] 
#>   [4] Inverse Function Theorem:    [PASSED] 
#>   [5] Link Derivatives:            [PASSED] 
#>   [6] Inverse Link Derivatives:    [PASSED] 
```

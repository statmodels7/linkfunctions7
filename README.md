---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->
[![R-CMD-check](https://github.com/statmodels7/linkfunctions7/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/statmodels7/linkfunctions7/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/statmodels7/linkfunctions7/graph/badge.svg)](https://app.codecov.io/gh/statmodels7/linkfunctions7)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->



# linkfunctions7

The goal of `{linkfunctions7}` is to provide a collection of mathematical link functions using the S7 object-oriented programming system in R, providing exact analytical derivatives up to the fourth order for both the forward and inverse link functions.

## Installation

You can install the development version of `{linkfunctions7}` from GitHub with:


``` r
# install.packages("pak")
pak::pak("statmodels7/linkfunctions7")
```

## Overview and Examples
The package revolves around a core `link` class. It includes specialized link functions and generic methods to evaluate the forward link $\eta = g(\theta)$, the inverse link $\theta = g^{-1}(\eta)$, and their analytical derivatives.

### The Softplus Link
The `softplus_link()` is designed to model strictly positive parameters, enforcing the domain constraint $\theta > 0$. It serves as a robust alternative to the standard `log_link()`. While a log link implies a global exponential relationship, the Softplus transformation smoothly transitions to a linear asymptote for large positive values of the linear predictor $\eta$.

``` r
softplus_link(a = 1)
#> S7 Link Object: softplus(a=1)
#>   - Parameter domain (theta): (0, Inf)
#>   - Link parameters: a = 1
```

You can visualize the behavior of any link object over its valid domain using the provided `plot()` method:


``` r
plot(softplus_link(a = 1))
```

<div class="figure">
<img src="man/figures/README-unnamed-chunk-4-1.png" alt="plot of chunk unnamed-chunk-4" width="100%" />
<p class="caption">plot of chunk unnamed-chunk-4</p>
</div>


### Bounded Links and Derivatives
The `bounded_link()` function generates a link object that restricts $\theta$ to a specified interval. When both lower and upper bounds are provided (a doubly bounded link), the function maps the interval $[\text{lwr}, \text{upr}]$ to the whole real line.

The package exposes generic functions to evaluate the exact analytical derivatives for these transformations. The following example demonstrates how to evaluate the inverse link, followed by the first, second, and fourth order derivatives of the forward link function $g(\theta)$.

``` r
link <- bounded_link(lwr = -3, upr = 2)
eta <- seq(-3, 3, l = 5)

theta <- linkinv(link, eta)
dlinkfun(link, theta)
#> [1] 4.427065 1.340964 0.800000 1.340964 4.427065
d2linkfun(link, theta)
#> [1] -17.739913  -1.142115   0.000000   1.142115  17.739913
d4linkfun(link, theta)
#> [1] -1897.611144    -8.646712     0.000000     8.646712  1897.611144
```

## Diagnostic Validation
To maintain mathematical rigor, `{linkfunctions7}` includes a comprehensive diagnostic method called `check_link()`. This function empirically verifies the algebraic invertibility of the transformations (confirming that $g^{-1}(g(\theta)) = \theta$). Furthermore, it tests strict monotonicity, the inverse function theorem, and validates all implemented analytical derivatives against numerical gradients to ensure accuracy.


``` r
check_link(log_link())
#> Checking S7 Link Object: log 
#>   [1] Invertibility (Theta space): [PASSED] 
#>   [2] Invertibility (Eta space):   [PASSED] 
#>   [3] Strict Monotonicity:         [PASSED] 
#>   [4] Inverse Function Theorem:    [PASSED] 
#>   [5] Link Derivatives:            [PASSED] 
#>   [6] Inverse Link Derivatives:    [PASSED]
```



<!-- README.md is generated from README.Rmd. Please edit that file, then
     regenerate with devtools::build_readme(). Do not use knitr::knit(): it
     processes the code but leaves this YAML header in the output as literal
     text, which GitHub and pkgdown both render verbatim. -->

<!-- badges: start -->

[![R-CMD-check](https://github.com/statmodels7/linkfunctions7/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/statmodels7/linkfunctions7/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/statmodels7/linkfunctions7/graph/badge.svg)](https://app.codecov.io/gh/statmodels7/linkfunctions7)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

# linkfunctions7 <img src="man/figures/logo.png" align="right" height="139" alt="" />

In most R modeling packages a link function is a character string passed
to a fitting routine, which turns it internally into a few closures.
Code outside that package cannot use those closures, compute their
higher derivatives, or extend them without editing the package source.

`{linkfunctions7}` makes a link an object. Each constructor in the
package returns a link with **exact analytical derivatives up to fifth
order in both directions**, forward and inverse, and a diagnostic
verifies those derivatives against numerical ones.

It is part of [statmodels7](https://statmodels7.github.io), an S7
toolkit for statistical modeling, and is what
[distributions7](https://statmodels7.github.io/distributions7/) uses to
move between a constrained parameter and the unconstrained scale a
fitting routine works on. The mathematics behind the formulas, and the
derivations the derivatives come from, is worked out in [the statmodels7
book](https://statmodels7.github.io/book/).

## Installation

``` r
# install.packages("pak")
pak::pak("statmodels7/linkfunctions7")
```

Or install the whole toolkit at once:

``` r
pak::pak("statmodels7/statmodels7")
```

## A link is an object

Each link is created by a constructor and carries its own name, domain
and parameters. The forward link $\eta = g(\theta)$ and the inverse
$\theta = g^{-1}(\eta)$ are generics that dispatch on it, as are all ten
derivatives.

The softplus link illustrates why more than one positivity link is
useful. Like the log it keeps $\theta$ positive, but where the log
implies an exponential relationship everywhere, the softplus bends
smoothly into a straight line as $\eta$ grows, so a large linear
predictor is not amplified into an enormous parameter:

``` r
softplus_link(a = 1)
#> S7 Link Object: softplus(a=1)
#>   - Parameter domain (theta): (0, Inf)
#>   - Link parameters: a = 1
plot(softplus_link(a = 1))
```

<img src="man/figures/README-unnamed-chunk-4-1.png" alt="" width="100%" />

A parameter confined to an interval gets `bounded_link()`, which maps
the interval onto the whole real line and accepts a lower bound, an
upper bound, or both. For every link in the package the derivatives are
closed-form expressions, not finite differences, and they go to fifth
order in both directions:

``` r
link <- bounded_link(lwr = -3, upr = 2)
theta <- linkinv(link, seq(-3, 3, length.out = 5))

theta                    # stays inside (-3, 2) whatever eta is
#> [1] -2.762871 -2.087872 -0.500000  1.087872  1.762871
dlinkfun(link, theta)    # g'(theta), exactly
#> [1] 4.427065 1.340964 0.800000 1.340964 4.427065
d5linkfun(link, theta)   # ... down to the fifth order
#> [1] 32009.92634    38.03412     0.49152    38.03412 32009.92634
```

## Validating a link

The derivatives are hand-written, so the package ships the tool that
checks them. `check_link()` confirms that a link inverts cleanly in both
directions, that it is strictly monotone, that the inverse function
theorem $g'(\theta)\,(g^{-1})'(\eta) = 1$ holds, and that every
analytical derivative agrees with a numerical one. It accepts any link
and is most useful for a newly written one:

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

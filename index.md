# linkfunctions7

In most R modelling packages a link function has no standing of its own.
It is a string passed to a fitting routine and unpacked internally into
a couple of closures that nothing outside can reach: it cannot be handed
to another package, asked for its second derivative, or extended without
editing the source that owns it.

[linkfunctions7](https://statmodels7.github.io/linkfunctions7/) makes a
link an object. The package provides fourteen constructors, each
carrying **exact analytical derivatives up to fourth order in both
directions**, forward and inverse, together with a diagnostic that
verifies those derivatives against numerical ones.

It is part of [statmodels7](https://statmodels7.github.io), an S7
toolkit for statistical modelling, and is what
[distributions7](https://statmodels7.github.io/distributions7) uses to
move between a constrained parameter and the unconstrained scale a
fitting routine works on. The mathematics behind every formula,
including the derivation of the derivatives to fourth order, is worked
out in [the statmodels7 book](https://statmodels7.github.io/book/).

## Installation

``` r

# install.packages("pak")
pak::pak("statmodels7/linkfunctions7")
```

## A link is an object

Each link is created by a constructor and knows its own name, domain and
parameters. The forward link $`\eta = g(\theta)`$ and the inverse
$`\theta = g^{-1}(\eta)`$ are generics that dispatch on it, as are all
eight derivatives.

The softplus link illustrates why more than one positivity link is
useful. Like the log it keeps $`\theta`$ positive, but where the log
implies an exponential relationship everywhere, the softplus bends
smoothly into a straight line as $`\eta`$ grows, so a large linear
predictor is not amplified into an enormous parameter:

``` r

softplus_link(a = 1)
#> S7 Link Object: softplus(a=1)
#>   - Parameter domain (theta): (0, Inf)
#>   - Link parameters: a = 1
plot(softplus_link(a = 1))
```

![](reference/figures/README-unnamed-chunk-3-1.png)

A parameter confined to an interval gets
[`bounded_link()`](https://statmodels7.github.io/linkfunctions7/reference/bounded_link.md),
which maps the interval onto the whole real line and accepts a lower
bound, an upper bound, or both. Whatever the link, the derivatives are
exact formulas rather than finite differences, and they go to fourth
order in both directions:

``` r

link <- bounded_link(lwr = -3, upr = 2)
theta <- linkinv(link, seq(-3, 3, length.out = 5))

theta                    # stays inside (-3, 2) whatever eta is
#> [1] -2.762871 -2.087872 -0.500000  1.087872  1.762871
dlinkfun(link, theta)    # g'(theta), exactly
#> [1] 4.427065 1.340964 0.800000 1.340964 4.427065
d4linkfun(link, theta)   # ... down to the fourth order
#> [1] -1897.611144    -8.646712     0.000000     8.646712  1897.611144
```

## Validating a link

The derivatives are hand-written, so the package ships the tool that
checks them.
[`check_link()`](https://statmodels7.github.io/linkfunctions7/reference/check_link.md)
confirms that a link inverts cleanly in both directions, that it is
strictly monotone, that the inverse function theorem
$`g'(\theta)\,(g^{-1})'(\eta) = 1`$ holds, and that every analytical
derivative agrees with a numerical one. It runs on any link, and it is
most useful on a newly written one, where a wrong derivative is
likeliest:

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

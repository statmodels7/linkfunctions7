# The Scalar Route of a Link

Returns how the package's compiled scalar entry points address a link:
the class name that the C function `lf7_class_id()` recognizes, and the
link's own parameters in the order `lf7_inv12p()` reads them. A consumer
that resolves the entry points with `R_GetCCallable()` uses the name to
obtain the link's identifier once and passes the parameters at every
evaluation.

## Usage

``` r
link_scalar_route(x)
```

## Arguments

- x:

  A link object inheriting from
  [`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md).

## Value

`NULL` when the entry points do not cover the link; otherwise a list
with `name`, a single character string, and `par`, a numeric vector of
the link's parameters (of length zero for a link without parameters).

## Details

The entry points are `lf7_class_id(name)` and
`lf7_inv12p(id, par, eta, h, h1, h2)`, which writes the inverse link and
its first two derivatives at one value of the linear predictor, each as
the class's R method computes it, and `lf7_clamp(theta, lwr, upr)`,
which applies
[`link_bounds_clamp()`](https://statmodels7.github.io/linkfunctions7/reference/link_bounds_clamp.md)
to one value. The parameters are \\\lambda\\ for
[`power_link()`](https://statmodels7.github.io/linkfunctions7/reference/power_link.md),
\\a\\ for
[`softplus_link()`](https://statmodels7.github.io/linkfunctions7/reference/softplus_link.md),
the lower bound and the width for a doubly bounded
[`bounded_link()`](https://statmodels7.github.io/linkfunctions7/reference/bounded_link.md),
and the single bound for a one-sided one; the other links have none.

The class is matched by its exact name, so a class defined outside the
package, including one that inherits from a class of the package and may
override its methods, returns `NULL`.

## Examples

``` r
link_scalar_route(logit_link())
#> $name
#> [1] "LogitLink"
#> 
#> $par
#> numeric(0)
#> 
link_scalar_route(bounded_link(lwr = 2, upr = 5))
#> $name
#> [1] "DoublyBoundedLink"
#> 
#> $par
#> [1] 2 3
#> 
link_scalar_route(power_link(0.5))
#> $name
#> [1] "PowerLink"
#> 
#> $par
#> [1] 0.5
#> 
```

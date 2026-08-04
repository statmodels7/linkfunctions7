# Is a Class the Base Link Class

Answers whether an S7 class is this package's own
[`link`](https://statmodels7.github.io/linkfunctions7/reference/link.md)
class, which is how a method inherited from the base class is told from
one a link registered for itself.

## Usage

``` r
is_base_link_class(cls)
```

## Arguments

- cls:

  An S7 class.

## Value

`TRUE` or `FALSE`.

## Details

Object identity is tried first, since it is the usual case and costs
nothing, and the class name and package are compared after. The second
comparison is what makes the answer reliable:
[`identical()`](https://rdrr.io/r/base/identical.html) on an S7 class is
object identity, so it is false for a class re-created from the same
definition, which is what happens whenever a package's code is evaluated
rather than loaded. A base fallback mistaken for an analytic method
makes every fallback differentiate the order below it, which is the
nested differencing the design exists to forbid.

## See also

[`link_fallback_orders`](https://statmodels7.github.io/linkfunctions7/reference/link_fallback_orders.md)

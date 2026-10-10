# Is a Class the Base Link Class

Returns whether an S7 class is this package's own
[`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md)
class. It is used to tell a method inherited from the base class from
one that a link registered for itself.

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
comparison is the one that makes the answer reliable.
[`identical()`](https://rdrr.io/r/base/identical.html) on an S7 class is
object identity, so it is false for a class re-created from the same
definition, which happens whenever a package's code is evaluated instead
of loaded, as it is under coverage instrumentation. If a base fallback
were mistaken for an analytic method, every fallback would differentiate
the order below it, which is the nested differencing that the fallbacks
are written to avoid.

## See also

[`link_fallback_orders()`](https://statmodels7.github.io/linkfunctions7/reference/link_fallback_orders.md)

# Which Derivative Orders a Link Computes Exactly

Reports, for each direction, how many derivative orders the link
implements analytically and which are therefore obtained by finite
differences.

## Usage

``` r
link_fallback_orders(x)
```

## Arguments

- x:

  An object of class `link`.

## Value

A list with `forward` and `inverse`, each an integer: the number of
leading orders implemented analytically, from 0 to 5.

## Details

Every link has a method for every derivative generic, because the base
class supplies numerical fallbacks for the orders that a link does not
implement. This function reports which orders are analytic. A fallback
is accurate only to the precision of a finite difference, so
[`check_link()`](https://statmodels7.github.io/linkfunctions7/reference/check_link.md)
reports such orders separately instead of counting them as passed.

## See also

[`check_link()`](https://statmodels7.github.io/linkfunctions7/reference/check_link.md)

## Examples

``` r
# everything the package ships is exact to fifth order
link_fallback_orders(logit_link())
#> $forward
#> [1] 5
#> 
#> $inverse
#> [1] 5
#> 
```

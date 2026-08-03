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
leading orders implemented analytically, from 0 to 4.

## Details

Every link can answer every derivative generic, because the base class
supplies numerical fallbacks for the orders a link does not implement.
That convenience requires a way of asking which is which — a fallback is
correct but not exact, and it is the reason
[`check_link`](https://statmodels7.github.io/linkfunctions7/reference/check_link.md)
reports such orders separately rather than passing them.

## See also

[`check_link`](https://statmodels7.github.io/linkfunctions7/reference/check_link.md)

## Examples

``` r
# everything the package ships is exact to fourth order
link_fallback_orders(logit_link())
#> $forward
#> [1] 4
#> 
#> $inverse
#> [1] 4
#> 
```

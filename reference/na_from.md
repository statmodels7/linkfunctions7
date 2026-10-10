# Carry Missingness From an Input Over to a Result

Sets `r` to `NA` wherever `v` is `NA`.

## Usage

``` r
na_from(r, v)
```

## Arguments

- r:

  A numeric vector, the computed result.

- v:

  The numeric vector the result was computed from.

## Value

`r`, with `NA` in every position where `v` is `NA`.

## Details

This handles the same problem as
[`const_like()`](https://statmodels7.github.io/linkfunctions7/reference/const_like.md)
for a computed result: an expression whose exponent vanishes no longer
depends on its argument, and loses the argument's missingness with it.
The power link is the case in point, since `theta^(lambda - 2)` is
exactly `1` for a missing `theta` once `lambda` is 2.

## See also

[`const_like()`](https://statmodels7.github.io/linkfunctions7/reference/const_like.md)

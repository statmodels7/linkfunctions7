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

Same hazard as
[`const_like`](https://statmodels7.github.io/linkfunctions7/reference/const_like.md),
one step further along: an expression whose exponent happens to vanish
stops depending on its argument, and loses the argument's missingness
along with it. The power link is where this bites, `theta^(lambda - 2)`
being exactly `1` for a missing `theta` once `lambda` is 2.

## See also

[`const_like`](https://statmodels7.github.io/linkfunctions7/reference/const_like.md)

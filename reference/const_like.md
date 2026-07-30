# A Constant Vector That Preserves Missingness

Returns `value` repeated to the length of `v`, but missing wherever `v`
is missing.

## Usage

``` r
const_like(v, value)
```

## Arguments

- v:

  A numeric vector whose length and missingness pattern are copied.

- value:

  The constant to repeat.

## Value

A numeric vector as long as `v`, equal to `value` except where `v` is
`NA`.

## Details

A derivative that reduces to a constant must still report that it does
not know the answer for an input it was not given. R makes this easy to
get wrong: `NA^0` is `1`, so `theta^(lambda - 2)` silently turns a
missing parameter into a number as soon as `lambda` is 2. Every
derivative method that returns a constant (the identity link's, the
square root link's third and fourth inverse derivatives) goes through
this helper instead of [`rep()`](https://rdrr.io/r/base/rep.html).

## See also

[`na_from`](https://statmodels7.github.io/linkfunctions7/reference/na_from.md),
the same idea for a computed result.

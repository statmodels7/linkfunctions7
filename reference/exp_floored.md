# A Floored Exponential

`exp(eta)`, bounded below by
[`exp_floor`](https://statmodels7.github.io/linkfunctions7/reference/exp_floor.md).

## Usage

``` r
exp_floored(eta)
```

## Arguments

- eta:

  A numeric vector of linear predictors.

## Value

A numeric vector, never smaller than
[`exp_floor`](https://statmodels7.github.io/linkfunctions7/reference/exp_floor.md).

## Details

A profile of a plain gaussian fit puts this `pmax` above the QR
decomposition of the same fit, which invites replacing it with a
[`min()`](https://rdrr.io/r/base/Extremes.html) reduction and an early
return. Measured, that is not worth doing: it gains 7 to 11 per cent on
the call itself and LOSES on the fit, because the reduction is a second
pass over the same vector and the allocation it avoids was not what the
profile was really charging for. The simple form is kept.

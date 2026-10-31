# Recompute a Power Term Without a Premature Overflow

Replaces the elements of `r`, a derivative of the form \\c\\v^q\\ that
was computed as a quotient or a product, by the same quantity formed
with the power split in two halves, \\(c\\v^{q/2})\\v^{q/2}\\, where `r`
is infinite or below 1e-290 (zero, or subnormal with its precision lost,
a loss that the coefficient then multiplies) although the value is not.

## Usage

``` r
rescue_pow(r, v, cf, q)
```

## Arguments

- r:

  A numeric vector, the first computation of the derivative.

- v:

  A numeric vector, the base of the power (\\\theta\\, \\\eta\\, or a
  distance to a bound).

- cf:

  The coefficient \\c\\.

- q:

  The exponent \\q\\ of the power.

## Value

A numeric vector of the same length as `r`.

## Details

A derivative such as \\-120/\theta^6\\ divides by a power that overflows
before the quotient underflows (\\\theta \> 10^{51.4}\\), and
\\24/\theta^5\\ multiplies a small coefficient by a power that overflows
before the product does. Only the elements that came out infinite or
tiny are recomputed, so the values of ordinary size are unchanged to the
last bit. An element that is still infinite or zero after the split is a
true overflow or underflow.

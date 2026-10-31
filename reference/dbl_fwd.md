# Forward Derivatives of the Doubly Bounded Link

Computes the derivative of order `k` of the doubly bounded link,
\\(k-1)!\\\[(-1)^{k-1}/d^k + 1/u^k\]\\, with \\d = \theta - l\\ and \\u
= W - d\\ the distances to the two bounds. It is the body of the six
forward derivative methods of the class.

## Usage

``` r
dbl_fwd(x, theta, k)
```

## Arguments

- x:

  A doubly bounded link.

- theta:

  A numeric vector inside the bounds.

- k:

  The derivative order, from 1 to 6.

## Value

A numeric vector of the same length as `theta`.

## Details

The distance \\d\\ is exact near the lower bound (Sterbenz) and \\u =
W - d\\ is exact near the upper one, where \\1 - p\\ with \\p = d/W\\
would cancel. The powers are of the reciprocals, split in two halves
around the factorial, so that neither \\W^k\\ nor \\d^k\\ underflows or
overflows before the value does.

# A Finite-Difference Step for a Given Order

The step for a central stencil of order `order`, scaled by the magnitude
of the evaluation point and, when bounds are supplied, shrunk so that
the whole stencil stays strictly inside them.

## Usage

``` r
fd_step(x, order, bounds = NULL)
```

## Arguments

- x:

  A numeric vector of evaluation points.

- order:

  The derivative order, 1 to 4.

- bounds:

  An optional length-2 numeric vector, the open interval `x` must stay
  inside.

## Value

A numeric vector of steps, the same length as `x`.

## Details

The unclamped step is \\\varepsilon^{1/(k+2)}\max(1, \lvert x\rvert)\\,
which balances truncation against rounding for a central difference of
order \\k\\: higher orders divide by a higher power of \\h\\, so they
need a larger one.

The clamp matters because a link's domain is open. Differentiating the
log link at \\\theta = 10^{-8}\\ with a step chosen from the magnitude
alone would evaluate \\\log\\ at a negative number; the stencil for
order 3 and 4 reaches \\2h\\, so it is that reach, not \\h\\, that must
fit.

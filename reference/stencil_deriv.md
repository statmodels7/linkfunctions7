# One Central Stencil, Never Nested

The order-`order` central finite difference of `f` at `x`, applied in a
single step rather than by composing lower-order differences.

## Usage

``` r
stencil_deriv(f, x, order, h)
```

## Arguments

- f:

  A vectorized function of one numeric argument.

- x:

  A numeric vector of evaluation points.

- order:

  The derivative order, 1 to 5.

- h:

  A numeric vector of steps, from
  [`fd_step()`](https://statmodels7.github.io/linkfunctions7/reference/fd_step.md).

## Value

A numeric vector of the same length as `x`.

## Details

The stencils of orders one to four are \$\$f' \approx \frac{f(x+h) -
f(x-h)}{2h}, \qquad f'' \approx \frac{f(x+h) - 2f(x) +
f(x-h)}{h^{2}},\$\$ \$\$f''' \approx \frac{f(x+2h) - 2f(x+h) + 2f(x-h) -
f(x-2h)}{2h^{3}}, \qquad f'''' \approx \frac{f(x+2h) - 4f(x+h) + 6f(x) -
4f(x-h) + f(x-2h)}{h^{4}}.\$\$

Applying one stencil of order \\k\\ differs from applying \\k\\ stencils
of order one, and this function exists because of the difference: each
numerical differentiation amplifies the error of the one before it, so a
fourth derivative computed by four nested first differences is dominated
by rounding error. For the identity link, whose third derivative is
exactly zero, nested differentiation returns non-zero values.

# One Central Stencil, Never Nested

The order-`order` central finite difference of `f` at `x`, applied in a
single step rather than by composing lower-order differences.

## Usage

``` r
stencil_deriv(f, x, order, h)
```

## Arguments

- f:

  A vectorised function of one numeric argument.

- x:

  A numeric vector of evaluation points.

- order:

  The derivative order, 1 to 4.

- h:

  A numeric vector of steps, from
  [`fd_step`](https://statmodels7.github.io/linkfunctions7/reference/fd_step.md).

## Value

A numeric vector of the same length as `x`.

## Details

The four stencils are \$\$f' \approx \frac{f(x+h) - f(x-h)}{2h}, \qquad
f'' \approx \frac{f(x+h) - 2f(x) + f(x-h)}{h^{2}},\$\$ \$\$f''' \approx
\frac{f(x+2h) - 2f(x+h) + 2f(x-h) - f(x-2h)}{2h^{3}}, \qquad f''''
\approx \frac{f(x+2h) - 4f(x+h) + 6f(x) - 4f(x-h) + f(x-2h)}{h^{4}}.\$\$

Applying one stencil of order \\k\\ is not the same as applying \\k\\
stencils of order one, and the difference is the whole reason this
function exists: each numerical differentiation multiplies the error of
the one before it, so a fourth derivative reached by four nested first
differences is noise. The identity link makes the point without any
arithmetic – its third derivative is exactly zero, and nested
differentiation returns a number of order one.

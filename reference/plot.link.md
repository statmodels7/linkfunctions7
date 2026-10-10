# Visualize Link Functions

Plot method for `link` objects. It generates a panel with two plots:

1.  The link function \\\eta = g(\theta)\\ over its valid domain.

2.  The inverse link function \\\theta = g^{-1}(\eta)\\ over a standard
    range of linear predictors.

## Usage

``` r
# S3 method for class 'link'
plot(x, ...)
```

## Arguments

- x:

  An object of class `link`.

- ...:

  Named graphical parameters passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html) for
  both panels, where they replace the defaults (for example `col`, `lwd`
  or `main`).

## Value

No return value, called for side effects (plotting).

## Details

The plotting ranges depend on whether the link bounds are finite. The
function sets the graphical parameters (`par`) for a side-by-side layout
and restores the original settings on exit.

## Examples

``` r
plot(logit_link())

plot(softplus_link(2))

```

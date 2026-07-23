# Visualize Link Functions

A robust S7 plot method for objects of class `link`. It generates a
panel with two plots:

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

  Additional graphical parameters passed to
  [`plot`](https://rdrr.io/r/graphics/plot.default.html).

## Value

No return value, called for side effects (plotting).

## Details

The function automatically determines sensible plotting ranges based on
whether the link bounds are finite or infinite. It temporarily modifies
the graphical parameters (`par`) to create a side-by-side layout and
restores the original settings upon exit.

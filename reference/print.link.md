# Print Method for S7 Link Objects

Prints the name of a link, the open interval in which its parameter
lies, and its link parameters, if it has any. The output has two lines,
and a third for a link with parameters.

## Usage

``` r
# S3 method for class 'link'
print(x, ...)
```

## Arguments

- x:

  An object of class `link`.

- ...:

  Additional arguments passed to methods, currently unused.

## Value

`x`, invisibly. Called for the printed output.

## Details

The domain is printed as an open interval, `(0, 1)` for a probability
link and never `[0, 1]`, because a link never returns an endpoint;
[`link_bounds_clamp()`](https://statmodels7.github.io/linkfunctions7/reference/link_bounds_clamp.md)
ensures this in double precision.

The parameter line appears only for a link that has parameters. It names
them, so `power(lambda=2)` and `bounded(lwr=0, upr=10)` show the values
passed to the constructor.

## Examples

``` r
print(logit_link())
#> S7 Link Object: logit
#>   - Parameter domain (theta): (0, 1)

# links carrying parameters report them too
print(power_link(2))
#> S7 Link Object: power(lambda=2)
#>   - Parameter domain (theta): (0, Inf)
#>   - Link parameters: lambda = 2
print(bounded_link(0, 10))
#> S7 Link Object: bounded(lwr=0, upr=10)
#>   - Parameter domain (theta): (0, 10)
#>   - Link parameters: lwr = 0, upr = 10
```

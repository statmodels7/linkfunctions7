# Print Method for S7 Link Objects

A standard S7 print method for objects of class `link`. It displays the
name of the link function, its valid parameter domain, and any
additional parameters it may have (e.g., lambda for a power link).

## Usage

``` r
# S3 method for class 'link'
print(x, ...)
```

## Arguments

- x:

  An object of class `link`.

- ...:

  Additional arguments passed to methods (currently unused).

## Value

The function returns `x` invisibly.

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

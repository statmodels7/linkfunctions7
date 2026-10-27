# linkfunctions7: An S7 Framework for Link Functions

Link functions as objects, built on the S7 object system. Each link
provides the forward and inverse maps and their exact analytical
derivatives up to the fifth order in both directions. A user-defined
link needs only the two maps, and the package computes the missing
derivatives numerically from the highest analytical derivative that the
link supplies. A validator checks invertibility, monotonicity, the
inverse function theorem and every analytical derivative against a
numerical one, and reports which orders are numerical. Part of the
'statmodels7' toolkit.

## See also

Useful links:

- <https://statmodels7.github.io/linkfunctions7/>

- <https://github.com/statmodels7/linkfunctions7>

- Report bugs at <https://github.com/statmodels7/linkfunctions7/issues>

## Author

**Maintainer**: Giovanni Tinervia <giovannitinervia9@gmail.com>

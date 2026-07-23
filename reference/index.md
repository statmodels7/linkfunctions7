# Package index

## The link class

The S7 class every link inherits from, and the diagnostic that checks
one.

- [`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md)
  : S7 Class for Statistical Link Functions
- [`check_link()`](https://statmodels7.github.io/linkfunctions7/reference/check_link.md)
  [`check_link.link()`](https://statmodels7.github.io/linkfunctions7/reference/check_link.md)
  : Validate and Check a Link Object

## Link functions

Each constructor returns a link object. Those taking arguments cover a
family: a power, a softplus scale, a bounded interval.

- [`identity_link()`](https://statmodels7.github.io/linkfunctions7/reference/identity_link.md)
  : The Identity Link Function
- [`log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.md)
  : The Logarithmic Link Function
- [`logit_link()`](https://statmodels7.github.io/linkfunctions7/reference/logit_link.md)
  : The Logit Link Function
- [`probit_link()`](https://statmodels7.github.io/linkfunctions7/reference/probit_link.md)
  : The Probit Link Function
- [`cloglog_link()`](https://statmodels7.github.io/linkfunctions7/reference/cloglog_link.md)
  : The Complementary Log-Log (ClogLog) Link Function
- [`loglog_link()`](https://statmodels7.github.io/linkfunctions7/reference/loglog_link.md)
  : The Log-Log Link Function
- [`cauchit_link()`](https://statmodels7.github.io/linkfunctions7/reference/cauchit_link.md)
  : The Cauchit Link Function
- [`rhobit_link()`](https://statmodels7.github.io/linkfunctions7/reference/rhobit_link.md)
  : The Rhobit (Fisher's z) Link Function
- [`sqrt_link()`](https://statmodels7.github.io/linkfunctions7/reference/sqrt_link.md)
  : The Square Root Link Function
- [`inverse_link()`](https://statmodels7.github.io/linkfunctions7/reference/inverse_link.md)
  : The Inverse (Reciprocal) Link Function
- [`inverse_sq_link()`](https://statmodels7.github.io/linkfunctions7/reference/inverse_sq_link.md)
  : The Inverse Square Link Function
- [`power_link()`](https://statmodels7.github.io/linkfunctions7/reference/power_link.md)
  : The Power Link Function
- [`softplus_link()`](https://statmodels7.github.io/linkfunctions7/reference/softplus_link.md)
  : The Softplus Link Function
- [`bounded_link()`](https://statmodels7.github.io/linkfunctions7/reference/bounded_link.md)
  : The General Bounded Link Function

## Evaluating a link

The two directions, and any derivative order of either.

- [`linkfun()`](https://statmodels7.github.io/linkfunctions7/reference/linkfun.md)
  : Evaluate Forward Link Function
- [`linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/linkinv.md)
  : Evaluate Inverse Link Function
- [`linkderiv()`](https://statmodels7.github.io/linkfunctions7/reference/linkderiv.md)
  [`linkderiv.link()`](https://statmodels7.github.io/linkfunctions7/reference/linkderiv.md)
  : Evaluate Derivative of Link Function by Order
- [`linkinvderiv()`](https://statmodels7.github.io/linkfunctions7/reference/linkinvderiv.md)
  [`linkinvderiv.link()`](https://statmodels7.github.io/linkfunctions7/reference/linkinvderiv.md)
  : Evaluate Derivative of Inverse Link Function by Order

## Derivatives by order

The order-specific generics.
[`linkderiv()`](https://statmodels7.github.io/linkfunctions7/reference/linkderiv.md)
and
[`linkinvderiv()`](https://statmodels7.github.io/linkfunctions7/reference/linkinvderiv.md)
route to these; call them directly where the extra dispatch matters.

- [`dlinkfun()`](https://statmodels7.github.io/linkfunctions7/reference/dlinkfun.md)
  : 1st Derivative of Link Function
- [`d2linkfun()`](https://statmodels7.github.io/linkfunctions7/reference/d2linkfun.md)
  : 2nd Derivative of Link Function
- [`d3linkfun()`](https://statmodels7.github.io/linkfunctions7/reference/d3linkfun.md)
  : 3rd Derivative of Link Function
- [`d4linkfun()`](https://statmodels7.github.io/linkfunctions7/reference/d4linkfun.md)
  : 4th Derivative of Link Function
- [`dlinkinv()`](https://statmodels7.github.io/linkfunctions7/reference/dlinkinv.md)
  : 1st Derivative of Inverse Link Function
- [`d2linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/d2linkinv.md)
  : 2nd Derivative of Inverse Link Function
- [`d3linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/d3linkinv.md)
  : 3rd Derivative of Inverse Link Function
- [`d4linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/d4linkinv.md)
  : 4th Derivative of Inverse Link Function

## Methods

- [`print(`*`<link>`*`)`](https://statmodels7.github.io/linkfunctions7/reference/print.link.md)
  : Print Method for S7 Link Objects
- [`plot(`*`<link>`*`)`](https://statmodels7.github.io/linkfunctions7/reference/plot.link.md)
  : Visualize Link Functions

# Package index

## The link class

The S7 class every link inherits from, the diagnostic that checks one,
and the report of which derivative orders it computes exactly.

- [`link()`](https://statmodels7.github.io/linkfunctions7/reference/link.md)
  : S7 Class for Statistical Link Functions
- [`check_link()`](https://statmodels7.github.io/linkfunctions7/reference/check_link.md)
  [`check_link.link()`](https://statmodels7.github.io/linkfunctions7/reference/check_link.md)
  : Validate and Check a Link Object
- [`link_fallback_orders()`](https://statmodels7.github.io/linkfunctions7/reference/link_fallback_orders.md)
  : Which Derivative Orders a Link Computes Exactly

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
[`linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/linkinv.md)
guarantees a result strictly inside the parameter’s domain, which is
what lets its output be handed straight back to
[`linkfun()`](https://statmodels7.github.io/linkfunctions7/reference/linkfun.md)
or to a density that validates against open intervals.

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
- [`link_bounds_clamp()`](https://statmodels7.github.io/linkfunctions7/reference/link_bounds_clamp.md)
  : Clamp a Parameter Strictly Inside Its Domain

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

- [`CauchitLink()`](https://statmodels7.github.io/linkfunctions7/reference/CauchitLink.md)
  : S7 Class for the Cauchit Link
- [`ClogLogLink()`](https://statmodels7.github.io/linkfunctions7/reference/ClogLogLink.md)
  : S7 Class for the Complementary Log-Log Link
- [`DoublyBoundedLink()`](https://statmodels7.github.io/linkfunctions7/reference/DoublyBoundedLink.md)
  : S7 Class for a Doubly Bounded Link
- [`IdentityLink()`](https://statmodels7.github.io/linkfunctions7/reference/IdentityLink.md)
  : S7 Class for the Identity Link
- [`InverseLink()`](https://statmodels7.github.io/linkfunctions7/reference/InverseLink.md)
  : S7 Class for the Inverse Link
- [`InverseSqLink()`](https://statmodels7.github.io/linkfunctions7/reference/InverseSqLink.md)
  : S7 Class for the InverseSq Link
- [`LogLink()`](https://statmodels7.github.io/linkfunctions7/reference/LogLink.md)
  : S7 Class for the Logarithmic Link
- [`LogLogLink()`](https://statmodels7.github.io/linkfunctions7/reference/LogLogLink.md)
  : S7 Class for the LogLog Link
- [`LogitLink()`](https://statmodels7.github.io/linkfunctions7/reference/LogitLink.md)
  : S7 Class for the Logit Link
- [`LowerBoundedLink()`](https://statmodels7.github.io/linkfunctions7/reference/LowerBoundedLink.md)
  : S7 Class for a Lower Bounded Link
- [`PowerLink()`](https://statmodels7.github.io/linkfunctions7/reference/PowerLink.md)
  : S7 Class for the Power Link
- [`ProbitLink()`](https://statmodels7.github.io/linkfunctions7/reference/ProbitLink.md)
  : S7 Class for the Probit Link
- [`RhobitLink()`](https://statmodels7.github.io/linkfunctions7/reference/RhobitLink.md)
  : S7 Class for the Rhobit Link
- [`SoftplusLink()`](https://statmodels7.github.io/linkfunctions7/reference/SoftplusLink.md)
  : S7 Class for the Softplus Link
- [`SqrtLink()`](https://statmodels7.github.io/linkfunctions7/reference/SqrtLink.md)
  : S7 Class for the Sqrt Link
- [`UpperBoundedLink()`](https://statmodels7.github.io/linkfunctions7/reference/UpperBoundedLink.md)
  : S7 Class for an Upper Bounded Link
- [`analytic_order()`](https://statmodels7.github.io/linkfunctions7/reference/analytic_order.md)
  : Highest Analytically Implemented Derivative Order
- [`check_link()`](https://statmodels7.github.io/linkfunctions7/reference/check_link.md)
  [`check_link.link()`](https://statmodels7.github.io/linkfunctions7/reference/check_link.md)
  : Validate and Check a Link Object
- [`const_like()`](https://statmodels7.github.io/linkfunctions7/reference/const_like.md)
  : A Constant Vector That Preserves Missingness
- [`eta_bounds()`](https://statmodels7.github.io/linkfunctions7/reference/eta_bounds.md)
  : The Range a Stencil May Evaluate the Inverse Link On
- [`exp_floor`](https://statmodels7.github.io/linkfunctions7/reference/exp_floor.md)
  : The Smallest Parameter Value the Exponential Links Will Report
- [`exp_floored()`](https://statmodels7.github.io/linkfunctions7/reference/exp_floored.md)
  : A Floored Exponential
- [`fallback_deriv()`](https://statmodels7.github.io/linkfunctions7/reference/fallback_deriv.md)
  : The Body Shared by Every Numerical Fallback
- [`fd_step()`](https://statmodels7.github.io/linkfunctions7/reference/fd_step.md)
  : A Finite-Difference Step for a Given Order
- [`is_base_link_class()`](https://statmodels7.github.io/linkfunctions7/reference/is_base_link_class.md)
  : Is a Class the Base Link Class
- [`linkderiv()`](https://statmodels7.github.io/linkfunctions7/reference/linkderiv.md)
  [`linkderiv.link()`](https://statmodels7.github.io/linkfunctions7/reference/linkderiv.md)
  : Evaluate Derivative of Link Function by Order
- [`linkinvderiv()`](https://statmodels7.github.io/linkfunctions7/reference/linkinvderiv.md)
  [`linkinvderiv.link()`](https://statmodels7.github.io/linkfunctions7/reference/linkinvderiv.md)
  : Evaluate Derivative of Inverse Link Function by Order
- [`logistic_deriv()`](https://statmodels7.github.io/linkfunctions7/reference/logistic_deriv.md)
  : Derivatives of the Standard Logistic Function
- [`na_from()`](https://statmodels7.github.io/linkfunctions7/reference/na_from.md)
  : Carry Missingness From an Input Over to a Result
- [`stencil_deriv()`](https://statmodels7.github.io/linkfunctions7/reference/stencil_deriv.md)
  : One Central Stencil, Never Nested
- [`print(`*`<link>`*`)`](https://statmodels7.github.io/linkfunctions7/reference/print.link.md)
  : Print Method for S7 Link Objects
- [`plot(`*`<link>`*`)`](https://statmodels7.github.io/linkfunctions7/reference/plot.link.md)
  : Visualize Link Functions

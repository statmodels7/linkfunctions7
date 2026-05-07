#' @title Evaluate Forward Link Function
#' @include link_class.R
#' @param x An object of class \code{link}.
#' @param theta A numeric vector of parameters.
#' @return A numeric vector of the linear predictor.
#' @export
linkfun <- S7::new_generic("linkfun", "x", fun = function(x, theta) S7::S7_dispatch())

#' @title Evaluate Inverse Link Function
#' @param x An object of class \code{link}.
#' @param eta A numeric vector of linear predictors.
#' @return A numeric vector of probabilities/means.
#' @export
linkinv <- S7::new_generic("linkinv", "x", fun = function(x, eta) S7::S7_dispatch())

#' @title 1st Derivative of Link Function
#' @param x An object of class \code{link}.
#' @param theta A numeric vector.
#' @export
dlinkfun <- S7::new_generic("dlinkfun", "x", fun = function(x, theta) S7::S7_dispatch())

#' @title 2nd Derivative of Link Function
#' @param x An object of class \code{link}.
#' @param theta A numeric vector.
#' @export
d2linkfun <- S7::new_generic("d2linkfun", "x", fun = function(x, theta) S7::S7_dispatch())

#' @title 3rd Derivative of Link Function
#' @param x An object of class \code{link}.
#' @param theta A numeric vector.
#' @export
d3linkfun <- S7::new_generic("d3linkfun", "x", fun = function(x, theta) S7::S7_dispatch())

#' @title 4th Derivative of Link Function
#' @param x An object of class \code{link}.
#' @param theta A numeric vector.
#' @export
d4linkfun <- S7::new_generic("d4linkfun", "x", fun = function(x, theta) S7::S7_dispatch())

#' @title 1st Derivative of Inverse Link Function
#' @param x An object of class \code{link}.
#' @param eta A numeric vector.
#' @export
dlinkinv <- S7::new_generic("dlinkinv", "x", fun = function(x, eta) S7::S7_dispatch())

#' @title 2nd Derivative of Inverse Link Function
#' @param x An object of class \code{link}.
#' @param eta A numeric vector.
#' @export
d2linkinv <- S7::new_generic("d2linkinv", "x", fun = function(x, eta) S7::S7_dispatch())

#' @title 3rd Derivative of Inverse Link Function
#' @param x An object of class \code{link}.
#' @param eta A numeric vector.
#' @export
d3linkinv <- S7::new_generic("d3linkinv", "x", fun = function(x, eta) S7::S7_dispatch())

#' @title 4th Derivative of Inverse Link Function
#' @param x An object of class \code{link}.
#' @param eta A numeric vector.
#' @export
d4linkinv <- S7::new_generic("d4linkinv", "x", fun = function(x, eta) S7::S7_dispatch())
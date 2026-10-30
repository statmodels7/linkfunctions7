#' @title Print Method for S7 Link Objects
#'
#' @include generics.R
#' @include link_class.R
#' @description
#' Prints the name of a link, the open interval in which its parameter lies,
#' and its link parameters, if it has any. The output has two lines, and a
#' third for a link with parameters.
#'
#' @details
#' The domain is printed as an open interval, `(0, 1)` for a probability
#' link and never `[0, 1]`, because a link never returns an endpoint;
#' [link_bounds_clamp()] ensures this in double precision.
#'
#' The parameter line appears only for a link that has parameters. It names
#' them, so `power(lambda=2)` and `bounded(lwr=0, upr=10)` show the values
#' passed to the constructor.
#'
#' @param x An object of class `link`.
#' @param ... Additional arguments passed to methods, currently unused.
#'
#' @return `x`, invisibly. The function is called for its printed output.
#'
#' @examples
#' print(logit_link())
#'
#' # links carrying parameters report them too
#' print(power_link(2))
#' print(bounded_link(0, 10))
#'
#' @rdname print.link
#' @usage \method{print}{link}(x, ...)
#' @aliases print.link
#' @export
print.link <- function(x, ...) {
  cat(
    "S7 Link Object: ", x@link_name, "\n",
    "  - Parameter domain (theta): (", x@link_bounds[1], ", ", x@link_bounds[2], ")\n",
    sep = ""
  )

  # Intelligently print additional parameters if they exist
  if (!is.null(x@link_params) && length(x@link_params) > 0) {
    params_str <- paste(names(x@link_params), x@link_params, sep = " = ", collapse = ", ")
    cat("  - Link parameters: ", params_str, "\n", sep = "")
  }

  invisible(x)
}

S7::method(print, link) <- print.link

#' @title Visualize Link Functions
#'
#' @description
#' `autoplot()` returns a graph of a `link` object in two panels, and `plot()`
#' prints that graph:
#'
#' 1. the link function \eqn{\eta = g(\theta)} over its valid domain;
#' 2. the inverse link function \eqn{\theta = g^{-1}(\eta)} over a standard
#'    range of linear predictors.
#'
#' Dotted grey lines mark zero on both axes of each panel.
#'
#' @details
#' The plotting ranges depend on whether the link bounds are finite. A
#' bounded domain is drawn up to a distance of `1e-4` from each finite bound,
#' and an unbounded side extends five units past the other bound (or from
#' \eqn{-5} to \eqn{5} when neither bound is finite). The inverse is drawn
#' over the values of \eqn{\eta} reached in the first panel, restricted to
#' \eqn{[-6, 6]} when the two ranges overlap.
#'
#' The graph is drawn with \pkg{ggplot2} and its two panels are joined with
#' \pkg{patchwork}; both packages are suggested, and the methods signal an
#' error that names them when they are not installed. The returned object
#' can be modified like any other, for instance with
#' `+ patchwork::plot_annotation(title = "my title")`.
#'
#' @param x,object An object of class `link`.
#' @param ... Unused; an argument passed here signals an error.
#'
#' @return `autoplot()` returns a patchwork object holding the two ggplot2
#'   panels; `plot()` prints the same object and returns it invisibly.
#'
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("patchwork", quietly = TRUE)
#' plot(logit_link())
#' plot(softplus_link(2))
#'
#' @rdname plot.link
#' @usage \method{plot}{link}(x, ...)
#' @aliases plot.link autoplot.link
#' @export
plot.link <- function(x, ...) {
  p <- autoplot.link(x, ...)
  print(p)
  invisible(p)
}

#' @rdname plot.link
autoplot.link <- function(object, ...) {
  if (...length()) {
    stop("The graph of a link takes no further arguments; modify the ",
         "returned object instead.", call. = FALSE)
  }
  need_ggplot(patchwork = TRUE)
  x <- object

  lb <- x@link_bounds[1]
  ub <- x@link_bounds[2]
  eps <- 1e-4 # A smaller epsilon for better boundary visualization

  # Robustly define the theta sequence based on bounds
  if (is.finite(lb) && is.finite(ub)) {
    theta_seq <- seq(lb + eps, ub - eps, length.out = 1001)
  } else if (is.finite(lb) && !is.finite(ub)) {
    theta_seq <- seq(lb + eps, lb + 5, length.out = 1001)
  } else if (!is.finite(lb) && is.finite(ub)) {
    theta_seq <- seq(ub - 5, ub - eps, length.out = 1001)
  } else { # Both infinite
    theta_seq <- seq(-5, 5, length.out = 1001)
  }

  eta_vals <- linkfun(x, theta_seq)

  # Determine a sensible range for eta based on the link function's output
  # to avoid NaNs for links with restricted eta domains (e.g., inverse_sq_link)
  valid_eta <- eta_vals[is.finite(eta_vals)]
  if (length(valid_eta) > 0) {
    eta_min <- min(valid_eta)
    eta_max <- max(valid_eta)

    # Try to constrain within a standard [-6, 6] range for visual consistency,
    # but respect the actual valid domain of eta
    eta_start <- max(eta_min, -6)
    eta_end <- min(eta_max, 6)

    # If the natural range does not overlap with [-6, 6], use the natural range
    if (eta_start >= eta_end) {
      eta_start <- eta_min
      eta_end <- eta_max
    }

    # Fallback for constant eta (highly unlikely for a valid link function)
    if (abs(eta_end - eta_start) < 1e-6) {
      eta_start <- eta_start - 1
      eta_end <- eta_end + 1
    }
  } else {
    eta_start <- -5
    eta_end <- 5
  }

  eta_seq <- seq(eta_start, eta_end, length.out = 1001)

  # Suppress warnings gracefully in case evaluation hits undefined boundary areas
  theta_vals <- suppressWarnings(linkinv(x, eta_seq))

  panel <- function(h, v, title, xlab, ylab) {
    d <- data.frame(h = h, v = v)
    d <- d[is.finite(d$h) & is.finite(d$v), , drop = FALSE]
    ggplot2::ggplot(d, ggplot2::aes(x = .data$h, y = .data$v)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey50") +
      ggplot2::geom_vline(xintercept = 0, colour = "grey50") +
      ggplot2::geom_line(colour = "#1f5fa8", linewidth = 0.9) +
      ggplot2::labs(title = title, x = xlab, y = ylab) +
      ggplot2::theme_bw()
  }
  .data <- ggplot2::.data
  p1 <- panel(theta_seq, eta_vals, "link function",
              expression(theta), expression(eta == g(theta)))
  p2 <- panel(eta_seq, theta_vals, "inverse link function",
              expression(eta), expression(theta == g^{-1} * (eta)))
  patchwork::wrap_plots(p1, p2, nrow = 1L) +
    patchwork::plot_annotation(title = paste("Link:", x@link_name))
}

S7::method(plot, link) <- plot.link
S7::method(autoplot, link) <- autoplot.link

#' @title Link Derivative Wrapper
#' @description Routes to the correct forward derivative generic based on order.
#' @param x An object of class `link`.
#' @param theta A numeric vector.
#' @param order An integer (0 to 5).
#' @return A numeric vector of the same length as `theta`.
#' @rdname linkderiv
#' @keywords internal
linkderiv.link <- function(x, theta, order = 1) {
  # switch() on the integer directly: converting it to a character first costs
  # more than the branch it selects, on a function called once per parameter per
  # order by anything that works on the link scale.
  if (length(order) != 1L || is.na(order) || order < 0 || order > 5) {
    stop("The forward derivative order is not supported: it must be one of 0 to 5.", call. = FALSE)
  }
  switch(as.integer(order) + 1L,
         linkfun(x, theta),
         dlinkfun(x, theta),
         d2linkfun(x, theta),
         d3linkfun(x, theta),
         d4linkfun(x, theta),
         d5linkfun(x, theta)
  )
}
S7::method(linkderiv, link) <- linkderiv.link

#' @title Inverse Link Derivative Wrapper
#' @description Routes to the correct inverse derivative generic based on order.
#' @param x An object of class `link`.
#' @param eta A numeric vector.
#' @param order An integer (0 to 5).
#' @return A numeric vector of the same length as `eta`.
#' @rdname linkinvderiv
#' @keywords internal
linkinvderiv.link <- function(x, eta, order = 1) {
  if (length(order) != 1L || is.na(order) || order < 0 || order > 5) {
    stop("The inverse derivative order is not supported: it must be one of 0 to 5.", call. = FALSE)
  }
  switch(as.integer(order) + 1L,
         linkinv(x, eta),
         dlinkinv(x, eta),
         d2linkinv(x, eta),
         d3linkinv(x, eta),
         d4linkinv(x, eta),
         d5linkinv(x, eta)
  )
}
S7::method(linkinvderiv, link) <- linkinvderiv.link

#' @title Inverse Link Derivatives from the Order-Specific Generics
#' @description Evaluates [dlinkinv()] and the higher generics up to `order`
#'   and returns them as a list. It is the method of [linkinv_derivs()] for
#'   the base class.
#' @param x An object of class `link`.
#' @param eta A numeric vector.
#' @param order A whole number from 1 to 5.
#' @return A list of `order` numeric vectors of the same length as `eta`.
#' @rdname linkinv_derivs.link
#' @keywords internal
linkinv_derivs.link <- function(x, eta, order) {
  switch(as.integer(order),
         list(dlinkinv(x, eta)),
         list(dlinkinv(x, eta), d2linkinv(x, eta)),
         list(dlinkinv(x, eta), d2linkinv(x, eta), d3linkinv(x, eta)),
         list(dlinkinv(x, eta), d2linkinv(x, eta), d3linkinv(x, eta),
              d4linkinv(x, eta)),
         list(dlinkinv(x, eta), d2linkinv(x, eta), d3linkinv(x, eta),
              d4linkinv(x, eta), d5linkinv(x, eta))
  )
}
S7::method(linkinv_derivs, link) <- linkinv_derivs.link

#' @title Validate and Check a Link Object
#'
#' @description
#' Validates a `link` object numerically: invertibility in both
#' directions on a grid, strict monotonicity, the inverse function theorem
#' \eqn{h'(\eta)\,g'(\theta) = 1}, and every analytic derivative against one
#' numerical differentiation of the analytic order below it.
#'
#' @param x An object of class `link`.
#' @param tolerance Numeric tolerance for floating-point comparisons.
#' @param ... Additional arguments passed to methods.
#'
#' @details
#' The method performs six checks:
#'
#' 1. Invertibility (\eqn{\theta} space): verifies
#'    \eqn{g^{-1}(g(\theta)) = \theta} on a grid of parameter values.
#' 2. Invertibility (\eqn{\eta} space): verifies
#'    \eqn{g(g^{-1}(\eta)) = \eta} on a grid spanning the linear predictors
#'    that the link produces from the parameter grid.
#' 3. Strict monotonicity: checks that \eqn{g'(\theta)} has the same sign
#'    at every point of the grid, so that the map is one-to-one.
#' 4. Inverse function theorem: verifies
#'    \eqn{g'(\theta) \cdot (g^{-1})'(\eta) = 1}.
#' 5. Link derivatives: compares each analytic derivative of
#'    \eqn{g(\theta)}, up to the fifth order, with a numerical one.
#' 6. Inverse link derivatives: compares each analytic derivative of
#'    \eqn{g^{-1}(\eta)}, up to the fifth order, with a numerical one.
#'
#' Each analytic derivative of order \eqn{k} is compared with one numerical
#' differentiation of the analytic derivative of order \eqn{k - 1}, so the
#' errors of successive numerical differentiations do not compound.
#'
#' @importFrom numDeriv grad
#' @return Invisibly, a named list of the check results: the four scalar logicals
#'   `invertibility_theta`, `invertibility_eta`, `monotonicity` and
#'   `inverse_theorem`, plus `link_derivatives` and
#'   `inverse_link_derivatives`, each a logical vector of length five named
#'   `order_1` to `order_5`. In those two, `NA` means that the order was not
#'   checked: it is supplied by a numerical fallback, so the value and the
#'   reference would come from the same computation and would agree even for
#'   a wrong link. The number of
#'   orders actually implemented is carried on the result as the attribute
#'   `"analytic_orders"`; see [link_fallback_orders()]. A
#'   derivative that raises an error still counts as `FALSE`. The function is
#'   called mainly for the summary that it prints to the console.
#'
#' @examples
#' # every link the package ships passes all six checks
#' check_link(logit_link())
#'
#' res <- check_link(power_link(2))
#' res$link_derivatives
#' res$inverse_theorem
#'
#' # a user-defined link should pass the same checks
#' all(unlist(check_link(bounded_link(0, 10))))
#'
#' @rdname check_link
#' @keywords internal
check_link.link <- function(x, tolerance = 1e-5, ...) {

  cat("Checking S7 Link Object:", x@link_name, "\n")

  # 1. Generate evaluation points strictly inside valid bounds for theta.
  #
  # The inset has to leave room for the numerical differentiation performed
  # below. numDeriv's Richardson stencil reaches roughly 8e-4 * |x| away from
  # each point, so a grid coming within 1e-3 of the boundary is differentiated
  # using values from outside the domain: those come back NaN and the check
  # reports a failure for derivatives that are in fact exact. It is therefore a
  # fraction of the span being sampled rather than a fixed absolute distance.
  lb <- x@link_bounds[1]
  ub <- x@link_bounds[2]
  span <- if (all(is.finite(c(lb, ub)))) ub - lb else 5
  eps <- 0.02 * span

  if (is.finite(lb) && is.finite(ub)) {
    theta_seq <- seq(lb + eps, ub - eps, length.out = 15)
  } else if (is.finite(lb) && !is.finite(ub)) {
    theta_seq <- seq(lb + eps, lb + 5, length.out = 15)
  } else if (!is.finite(lb) && is.finite(ub)) {
    theta_seq <- seq(ub - 5, ub - eps, length.out = 15)
  } else {
    theta_seq <- seq(-3, 3, length.out = 15)
  }

  # 2. Test Algebraic Invertibility (Theta -> Eta -> Theta)
  eta_vals <- linkfun(x, theta_seq)
  theta_hat <- linkinv(x, eta_vals)

  inv_error <- max(abs(theta_seq - theta_hat))
  invertibility_pass <- !is.na(inv_error) && inv_error <= tolerance

  # 3. Test Algebraic Invertibility (Eta -> Theta -> Eta)
  #
  # Over the eta the link can actually produce, not a fixed [-4, 4]. A link whose
  # range is restricted -- the square root maps onto the positive half-line, the
  # inverse square onto (0, Inf) -- can say nothing sensible about an eta it never
  # produces, and testing there reported a failure for links that are perfectly
  # invertible on their own range.
  eta_seq_test <- seq(min(eta_vals), max(eta_vals), length.out = 15)
  eta_hat <- linkfun(x, linkinv(x, eta_seq_test))
  inv_eta_error <- max(abs(eta_seq_test - eta_hat))
  invertibility_eta_pass <- !is.na(inv_eta_error) && inv_eta_error <= tolerance

  # 4. Test Strict Monotonicity
  d1_theta <- linkderiv(x, theta_seq, order = 1)
  monotonicity_pass <- all(d1_theta > 0) || all(d1_theta < 0)

  # 5. Test Inverse Function Theorem (Derivative Reciprocal Identity)
  d1_eta <- linkinvderiv(x, eta_vals, order = 1)
  inv_thm_error <- max(abs(d1_eta * d1_theta - 1))
  inv_thm_pass <- !is.na(inv_thm_error) && inv_thm_error <= tolerance

  # Generic function to test a chain of derivatives
  test_derivative_chain <- function(eval_seq, deriv_fn, max_order = 5,
                                    n_exact = max_order) {
    results <- rep(NA, max_order)
    names(results) <- paste0("order_", 1:max_order)
    pass_prev <- TRUE

    for (o in 1:max_order) {
      if (!pass_prev) {
        results[o:max_order] <- FALSE
        break
      }

      # Beyond the orders the link implements itself there is nothing left to
      # check. The value would come from the numerical fallback, and the
      # reference is a numerical differentiation of the order below, so the two
      # are the same arithmetic and would agree however wrong the link is. That
      # is a vacuous pass, and a vacuous pass reported as [PASSED] is worse than
      # no check at all, so these orders are reported as numerical instead.
      if (o > n_exact) {
        results[o:max_order] <- NA
        break
      }

      fn_prev <- function(v) deriv_fn(x, v, order = o - 1)
      fn_curr <- function(v) deriv_fn(x, v, order = o)

      # A derivative that errors is a failure, not an absence. It used to leave
      # NA here and break out of the loop, and since the summary below reduces
      # with na.rm = TRUE, a link supplying only its first derivative reported
      # [PASSED] on all four orders -- the one result a user writing their own
      # link most needs to be told is wrong.
      curr_val <- tryCatch(fn_curr(eval_seq), error = function(e) NULL)
      if (is.null(curr_val)) {
        results[o:max_order] <- FALSE
        break
      }

      # numDeriv::grad() is element-wise when func maps a vector to a vector of
      # the same length, which is what a link's derivative does. One call
      # therefore replaces one call per evaluation point, for the same answer to
      # the last bit and about fifteen times faster.
      num_deriv <- numDeriv::grad(func = fn_prev, x = eval_seq)

      # Scale error relatively for large derivatives to avoid boundary inflation
      deriv_error <- max(abs(num_deriv - curr_val) / pmax(1, abs(curr_val)), na.rm = TRUE)
      if (is.na(deriv_error) || deriv_error > tolerance) {
        results[o] <- FALSE
        pass_prev <- FALSE
      } else {
        results[o] <- TRUE
      }
    }
    return(results)
  }

  # 6. Test Link Derivatives (Chained)
  exact <- link_fallback_orders(x)
  link_deriv_pass <- test_derivative_chain(theta_seq, linkderiv, max_order = 5,
                                           n_exact = exact$forward)

  # 7. Test Inverse Link Derivatives (Chained)
  inv_deriv_pass <- test_derivative_chain(eta_vals, linkinvderiv, max_order = 5,
                                          n_exact = exact$inverse)

  # A verdict for a family of orders: FAILED if any checked order failed,
  # PASSED if all checked orders passed, and neither if none was checked.
  deriv_verdict <- function(res, n_exact) {
    if (any(res %in% FALSE)) return("[FAILED]")
    if (n_exact == 0L) return("[numerical]")
    if (n_exact < length(res)) {
      return(sprintf("[PASSED to order %d, %d numerical]",
                     n_exact, length(res) - n_exact))
    }
    "[PASSED]"
  }

  # Concise console summary
  cat("  [1] Invertibility (Theta space):", if (invertibility_pass) "[PASSED]" else "[FAILED]", "\n")
  cat("  [2] Invertibility (Eta space):  ", if (invertibility_eta_pass) "[PASSED]" else "[FAILED]", "\n")
  cat("  [3] Strict Monotonicity:        ", if (monotonicity_pass) "[PASSED]" else "[FAILED]", "\n")
  cat("  [4] Inverse Function Theorem:   ", if (inv_thm_pass) "[PASSED]" else "[FAILED]", "\n")
  cat("  [5] Link Derivatives:           ", deriv_verdict(link_deriv_pass, exact$forward), "\n")
  cat("  [6] Inverse Link Derivatives:   ", deriv_verdict(inv_deriv_pass, exact$inverse), "\n")

  out <- list(
    invertibility_theta = invertibility_pass,
    invertibility_eta = invertibility_eta_pass,
    monotonicity = monotonicity_pass,
    inverse_theorem = inv_thm_pass,
    link_derivatives = link_deriv_pass,
    inverse_link_derivatives = inv_deriv_pass
  )
  # Carried as an attribute rather than a seventh element, so that every element
  # of the result stays a logical vector and callers can reduce over it.
  attr(out, "analytic_orders") <- exact
  invisible(out)
}
S7::method(check_link, link) <- check_link.link

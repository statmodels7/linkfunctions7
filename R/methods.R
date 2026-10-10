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
#' @return `x`, invisibly. Called for the printed output.
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
#' Plot method for `link` objects.
#' It generates a panel with two plots:
#'
#' 1. The link function \eqn{\eta = g(\theta)} over its valid domain.
#' 2. The inverse link function \eqn{\theta = g^{-1}(\eta)} over a standard range of linear predictors.
#'
#' @param x An object of class `link`.
#' @param ... Named graphical parameters passed to [graphics::plot()] for
#'   both panels, where they replace the defaults (for example `col`, `lwd`
#'   or `main`).
#'
#' @details
#' The plotting ranges depend on whether the link bounds are finite. The
#' function sets the graphical parameters (`par`) for a side-by-side layout
#' and restores the original settings on exit.
#'
#' @importFrom graphics par plot grid abline mtext
#'
#' @return No return value, called for side effects (plotting).
#'
#' @examples
#' plot(logit_link())
#' plot(softplus_link(2))
#'
#' @rdname plot.link
#' @usage \method{plot}{link}(x, ...)
#' @aliases plot.link
#' @export
plot.link <- function(x, ...) {
  dots <- list(...)
  if (length(dots) && (is.null(names(dots)) || any(names(dots) == ""))) {
    stop("Arguments passed through '...' must be named.", call. = FALSE)
  }
  # The caller's graphical parameters replace the defaults of each panel.
  draw <- function(defaults) {
    defaults[names(dots)] <- dots
    do.call(graphics::plot, defaults)
  }

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))

  graphics::par(
    mfrow = c(1, 2),
    mar = c(5, 5, 3.5, 1) + 0.1,
    oma = c(0, 0, 2.5, 0)
  )

  # --- Plot 1: Link function (eta vs. theta) ---

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

  draw(list(
    x = theta_seq, y = eta_vals,
    type = "l", lwd = 2, las = 1,
    xlab = expression(theta),
    ylab = expression(eta == g(theta)),
    main = "Link Function"
  ))
  graphics::grid()
  graphics::abline(h = 0, v = 0, lty = 3, col = "darkgray")

  # --- Plot 2: Inverse link function (theta vs. eta) ---

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

  draw(list(
    x = eta_seq, y = theta_vals,
    type = "l", lwd = 2, las = 1,
    xlab = expression(eta),
    ylab = expression(theta == g^{-1} * (eta)),
    main = "Inverse Link Function"
  ))
  graphics::grid()
  graphics::abline(h = 0, v = 0, lty = 3, col = "darkgray")

  graphics::mtext(
    text = paste("Link:", x@link_name),
    side = 3,
    line = 0.5,
    outer = TRUE,
    cex = 1.3,
    font = 2
  )
}

S7::method(plot, link) <- plot.link

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
    stop("Forward derivative order not supported.", call. = FALSE)
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
    stop("Inverse derivative order not supported.", call. = FALSE)
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
#' 1. **Invertibility (\eqn{\theta} space):** verifies
#'    \eqn{g^{-1}(g(\theta)) = \theta} on a grid of parameter values.
#' 2. **Invertibility (\eqn{\eta} space):** verifies
#'    \eqn{g(g^{-1}(\eta)) = \eta} on a grid spanning the linear predictors
#'    that the link produces from the parameter grid.
#' 3. **Strict monotonicity:** checks that \eqn{g'(\theta)} has the same sign
#'    at every point of the grid, so that the map is one-to-one.
#' 4. **Inverse function theorem:** verifies
#'    \eqn{g'(\theta) \cdot (g^{-1})'(\eta) = 1}.
#' 5. **Link derivatives:** compares each analytic derivative of
#'    \eqn{g(\theta)}, up to the fifth order, with a numerical one.
#' 6. **Inverse link derivatives:** compares each analytic derivative of
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
#'   derivative that raises an error still counts as `FALSE`. Called mainly
#'   for the summary printed to the console.
#'
#' @examples
#' # every link the package ships passes all six checks
#' check_link(logit_link())
#'
#' res <- check_link(power_link(2))
#' res$link_derivatives
#' res$inverse_theorem
#'
#' # the checks are what a user-defined link should be held to as well
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

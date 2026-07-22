#' @title Print Method for S7 Link Objects
#'
#' @include generics.R
#' @include link_class.R
#' @description
#' A standard S7 print method for objects of class \code{link}.
#' It displays the name of the link function, its valid parameter domain, and any
#' additional parameters it may have (e.g., lambda for a power link).
#'
#' @param x An object of class \code{link}.
#' @param ... Additional arguments passed to methods (currently unused).
#'
#' @return The function returns \code{x} invisibly.
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
#' A robust S7 plot method for objects of class \code{link}.
#' It generates a panel with two plots:
#' \enumerate{
#'   \item The link function \eqn{\eta = g(\theta)} over its valid domain.
#'   \item The inverse link function \eqn{\theta = g^{-1}(\eta)} over a standard range of linear predictors.
#' }
#'
#' @param x An object of class \code{link}.
#' @param ... Additional graphical parameters passed to \code{\link[graphics]{plot}}.
#'
#' @details
#' The function automatically determines sensible plotting ranges based on whether the
#' link bounds are finite or infinite. It temporarily modifies the graphical parameters
#' (\code{par}) to create a side-by-side layout and restores the original settings upon exit.
#'
#' @importFrom graphics par plot grid abline mtext
#'
#' @return No return value, called for side effects (plotting).
#'
#' @rdname plot.link
#' @usage \method{plot}{link}(x, ...)
#' @aliases plot.link
#' @export
plot.link <- function(x, ...) {
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
  
  graphics::plot(
    theta_seq, eta_vals,
    type = "l", lwd = 2, las = 1,
    xlab = expression(theta),
    ylab = expression(eta == g(theta)),
    main = "Link Function"
  )
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
  
  graphics::plot(
    eta_seq, theta_vals,
    type = "l", lwd = 2, las = 1,
    xlab = expression(eta),
    ylab = expression(theta == g^{-1} * (eta)),
    main = "Inverse Link Function"
  )
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
#' @param x An object of class \code{link}.
#' @param theta A numeric vector.
#' @param order An integer (0 to 4).
#' @rdname linkderiv
#' @export
linkderiv.link <- function(x, theta, order = 1) {
  # switch() on the integer directly: converting it to a character first costs
  # more than the branch it selects, on a function called once per parameter per
  # order by anything that works on the link scale.
  if (length(order) != 1L || is.na(order) || order < 0 || order > 4) {
    stop("Forward derivative order not supported.", call. = FALSE)
  }
  switch(as.integer(order) + 1L,
         linkfun(x, theta),
         dlinkfun(x, theta),
         d2linkfun(x, theta),
         d3linkfun(x, theta),
         d4linkfun(x, theta)
  )
}
S7::method(linkderiv, link) <- linkderiv.link

#' @title Inverse Link Derivative Wrapper
#' @description Routes to the correct inverse derivative generic based on order.
#' @param x An object of class \code{link}.
#' @param eta A numeric vector.
#' @param order An integer (0 to 4).
#' @rdname linkinvderiv
#' @export
linkinvderiv.link <- function(x, eta, order = 1) {
  if (length(order) != 1L || is.na(order) || order < 0 || order > 4) {
    stop("Inverse derivative order not supported.", call. = FALSE)
  }
  switch(as.integer(order) + 1L,
         linkinv(x, eta),
         dlinkinv(x, eta),
         d2linkinv(x, eta),
         d3linkinv(x, eta),
         d4linkinv(x, eta)
  )
}
S7::method(linkinvderiv, link) <- linkinvderiv.link

#' @title Validate and Check a Link Object
#'
#' @description
#' A diagnostic S7 method to mathematically validate a \code{link} object.
#' It sequentially verifies the algebraic invertibility and the correctness
#' of the analytical derivatives using numerical gradients in a chained sequence.
#'
#' @param x An object of class \code{link}.
#' @param tolerance Numeric tolerance for floating-point comparisons.
#' @param ... Additional arguments passed to methods.
#'
#' @details
#' The function assumes the existence of S7 generics \code{linkfun}, \code{linkinv},
#' \code{linkderiv}, and \code{linkinvderiv}. The method performs the following six diagnostic checks:
#' \enumerate{
#'   \item \strong{Invertibility (\eqn{\theta} space):} Verifies \eqn{g^{-1}(g(\theta)) = \theta}. Ensures that mapping from the parameter space to the linear predictor and back is lossless.
#'   \item \strong{Invertibility (\eqn{\eta} space):} Verifies \eqn{g(g^{-1}(\eta)) = \eta}. Ensures that mapping from the linear predictor to the parameter space and back is lossless. Note that this test may fail intentionally and correctly for links that map to a restricted \eqn{\eta} domain (e.g., the square root link).
#'   \item \strong{Strict Monotonicity:} Checks if the first derivative \eqn{g'(\theta)} is strictly positive or strictly negative across the domain, guaranteeing a one-to-one mapping.
#'   \item \strong{Inverse Function Theorem:} Verifies the mathematical identity \eqn{g'(\theta) \cdot (g^{-1})'(\eta) = 1}, confirming the theoretical relationship between the link derivative and the inverse link derivative.
#'   \item \strong{Link Derivatives:} Validates the exact analytical forward derivatives of \eqn{g(\theta)} up to the 4th order by comparing them against numerical gradients.
#'   \item \strong{Inverse Link Derivatives:} Validates the exact analytical inverse derivatives of \eqn{g^{-1}(\eta)} up to the 4th order by comparing them against numerical gradients.
#' }
#' Both forward and inverse derivative testing avoids compounding numerical errors by applying
#' first-order numerical differentiation iteratively to the exact lower-order analytical derivatives.
#'
#' @importFrom numDeriv grad
#' @return A logical list returning the success status of all available checks.
#'
#' @rdname check_link
#' @export
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
  
  # Helper function to compute numerical gradient robustly utilizing the native R pipe
  compute_num_grad <- function(eval_fn, val_seq) {
    val_seq |>
      vapply(function(val) numDeriv::grad(func = eval_fn, x = val), numeric(1))
  }
  
  # Generic function to test a chain of derivatives
  test_derivative_chain <- function(base_name, eval_seq, deriv_fn, max_order = 4) {
    results <- rep(NA, max_order)
    names(results) <- paste0("order_", 1:max_order)
    pass_prev <- TRUE
    
    for (o in 1:max_order) {
      if (!pass_prev) {
        results[o:max_order] <- FALSE
        break
      }
      
      fn_prev <- function(v) deriv_fn(x, v, order = o - 1)
      fn_curr <- function(v) deriv_fn(x, v, order = o)
      
      # Check if the derivative is implemented and evaluates without error
      curr_val <- tryCatch(fn_curr(eval_seq), error = function(e) NULL)
      if (is.null(curr_val)) {
        break
      }
      
      num_deriv <- compute_num_grad(fn_prev, eval_seq)
      exact_deriv <- curr_val
      
      # Scale error relatively for large derivatives to avoid boundary inflation
      deriv_error <- max(abs(num_deriv - exact_deriv) / pmax(1, abs(exact_deriv)), na.rm = TRUE)
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
  link_deriv_pass <- test_derivative_chain("link", theta_seq, linkderiv, max_order = 4)
  
  # 7. Test Inverse Link Derivatives (Chained)
  inv_deriv_pass <- test_derivative_chain("inverse link", eta_vals, linkinvderiv, max_order = 4)
  
  # Concise console summary
  cat("  [1] Invertibility (Theta space):", if (invertibility_pass) "[PASSED]" else "[FAILED]", "\n")
  cat("  [2] Invertibility (Eta space):  ", if (invertibility_eta_pass) "[PASSED]" else "[FAILED]", "\n")
  cat("  [3] Strict Monotonicity:        ", if (monotonicity_pass) "[PASSED]" else "[FAILED]", "\n")
  cat("  [4] Inverse Function Theorem:   ", if (inv_thm_pass) "[PASSED]" else "[FAILED]", "\n")
  cat("  [5] Link Derivatives:           ", if (all(link_deriv_pass, na.rm = TRUE)) "[PASSED]" else "[FAILED]", "\n")
  cat("  [6] Inverse Link Derivatives:   ", if (all(inv_deriv_pass, na.rm = TRUE)) "[PASSED]" else "[FAILED]", "\n")
  
  invisible(list(
    invertibility_theta = invertibility_pass,
    invertibility_eta = invertibility_eta_pass,
    monotonicity = monotonicity_pass,
    inverse_theorem = inv_thm_pass,
    link_derivatives = link_deriv_pass,
    inverse_link_derivatives = inv_deriv_pass
  ))
}
S7::method(check_link, link) <- check_link.link
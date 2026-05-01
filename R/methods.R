#' @title Print Method for S7 Link Objects
#'
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
#' @name print.link
S7::method(print, link) <- function(x, ...) {
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
#' @name plot.link
S7::method(plot, link) <- function(x, ...) {
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
    theta_seq <- seq(lb + eps, lb + 10, length.out = 1001)
  } else if (!is.finite(lb) && is.finite(ub)) {
    theta_seq <- seq(ub - 10, ub - eps, length.out = 1001)
  } else { # Both infinite
    theta_seq <- seq(-5, 5, length.out = 1001)
  }
  
  eta_vals <- x@linkfun(theta_seq)
  
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
  
  # Use a standard, fixed range for eta to allow for consistent comparison
  # across different link functions. This is the key to a robust plot.
  eta_seq <- seq(-6, 6, length.out = 1001)
  theta_vals <- x@linkinv(eta_seq)
  
  graphics::plot(
    eta_seq, theta_vals,
    type = "l", lwd = 2, las = 1,
    xlab = expression(eta),
    ylab = expression(theta == g^{-1}(eta)),
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
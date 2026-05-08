test_that("check_link mathematically validates implemented links", {
  # Capture the console output to keep the test summary clean
  # and suppress warnings that can arise at domain boundaries during checks.
  suppressWarnings(capture.output({
    res_bounded_double <- check_link(bounded_link(lwr = 0, upr = 1))
    res_bounded_lower <- check_link(bounded_link(lwr = 0))
    res_bounded_upper <- check_link(bounded_link(upr = 1))
    res_bounded_none <- check_link(bounded_link())
    res_cauchit <- check_link(cauchit_link())
    res_cloglog <- check_link(cloglog_link())
    res_identity <- check_link(identity_link())
    res_inverse <- check_link(inverse_link())
    res_inverse_sq <- check_link(inverse_sq_link())
    res_log <- check_link(log_link())
    res_logit <- check_link(logit_link())
    res_loglog <- check_link(loglog_link())
    res_power <- check_link(power_link())
    res_probit <- check_link(probit_link())
    res_rhobit <- check_link(rhobit_link())
    res_softplus <- check_link(softplus_link())
    res_sqrt <- check_link(sqrt_link())
  }))

  # Mathematical validations for bounded_link (double)
  expect_true(res_bounded_double$invertibility_theta, info = "Invertibility (theta) for bounded_double")
  expect_true(res_bounded_double$invertibility_eta, info = "Invertibility (eta) for bounded_double")
  expect_true(res_bounded_double$monotonicity, info = "Monotonicity for bounded_double")
  expect_true(res_bounded_double$inverse_theorem, info = "Inverse theorem for bounded_double")
  expect_true(all(res_bounded_double$link_derivatives, na.rm = TRUE), info = "Link derivatives for bounded_double")
  expect_true(all(res_bounded_double$inverse_link_derivatives, na.rm = TRUE), info = "Inverse link derivatives for bounded_double")

  # Mathematical validations for bounded_link (lower)
  expect_true(res_bounded_lower$invertibility_theta, info = "Invertibility (theta) for bounded_lower")
  expect_true(res_bounded_lower$invertibility_eta, info = "Invertibility (eta) for bounded_lower")
  expect_true(res_bounded_lower$monotonicity, info = "Monotonicity for bounded_lower")
  expect_true(res_bounded_lower$inverse_theorem, info = "Inverse theorem for bounded_lower")
  expect_true(all(res_bounded_lower$link_derivatives, na.rm = TRUE), info = "Link derivatives for bounded_lower")
  expect_true(all(res_bounded_lower$inverse_link_derivatives, na.rm = TRUE), info = "Inverse link derivatives for bounded_lower")

  # Mathematical validations for bounded_link (upper)
  expect_true(res_bounded_upper$invertibility_theta, info = "Invertibility (theta) for bounded_upper")
  expect_true(res_bounded_upper$invertibility_eta, info = "Invertibility (eta) for bounded_upper")
  expect_true(res_bounded_upper$monotonicity, info = "Monotonicity for bounded_upper")
  expect_true(res_bounded_upper$inverse_theorem, info = "Inverse theorem for bounded_upper")
  expect_true(all(res_bounded_upper$link_derivatives, na.rm = TRUE), info = "Link derivatives for bounded_upper")
  expect_true(all(res_bounded_upper$inverse_link_derivatives, na.rm = TRUE), info = "Inverse link derivatives for bounded_upper")

  # Mathematical validations for bounded_link (none)
  expect_true(res_bounded_none$invertibility_theta, info = "Invertibility (theta) for bounded_none")
  expect_true(res_bounded_none$invertibility_eta, info = "Invertibility (eta) for bounded_none")
  expect_true(res_bounded_none$monotonicity, info = "Monotonicity for bounded_none")
  expect_true(res_bounded_none$inverse_theorem, info = "Inverse theorem for bounded_none")
  expect_true(all(res_bounded_none$link_derivatives, na.rm = TRUE), info = "Link derivatives for bounded_none")
  expect_true(all(res_bounded_none$inverse_link_derivatives, na.rm = TRUE), info = "Inverse link derivatives for bounded_none")

  # Mathematical validations for cauchit_link
  expect_true(res_cauchit$invertibility_theta, info = "Invertibility (theta) for cauchit")
  expect_true(res_cauchit$invertibility_eta, info = "Invertibility (eta) for cauchit")
  expect_true(res_cauchit$monotonicity, info = "Monotonicity for cauchit")
  expect_true(res_cauchit$inverse_theorem, info = "Inverse theorem for cauchit")
  expect_true(all(res_cauchit$link_derivatives, na.rm = TRUE), info = "Link derivatives for cauchit")
  expect_true(all(res_cauchit$inverse_link_derivatives, na.rm = TRUE), info = "Inverse link derivatives for cauchit")

  # Mathematical validations for cloglog_link
  expect_true(res_cloglog$invertibility_theta, info = "Invertibility (theta) for cloglog")
  # Fails eta invertibility due to machine precision loss: 1 - exp(-exp(eta)) evaluates to exactly 1 for large eta.
  expect_false(res_cloglog$invertibility_eta, info = "Invertibility (eta) for cloglog (numerical precision)")
  expect_true(res_cloglog$monotonicity, info = "Monotonicity for cloglog")
  expect_true(res_cloglog$inverse_theorem, info = "Inverse theorem for cloglog")
  expect_true(all(res_cloglog$link_derivatives, na.rm = TRUE), info = "Link derivatives for cloglog")
  expect_true(all(res_cloglog$inverse_link_derivatives, na.rm = TRUE), info = "Inverse link derivatives for cloglog")

  # Mathematical validations for identity_link
  expect_true(res_identity$invertibility_theta, info = "Invertibility (theta) for identity")
  expect_true(res_identity$invertibility_eta, info = "Invertibility (eta) for identity")
  expect_true(res_identity$monotonicity, info = "Monotonicity for identity")
  expect_true(res_identity$inverse_theorem, info = "Inverse theorem for identity")
  expect_true(all(res_identity$link_derivatives, na.rm = TRUE), info = "Link derivatives for identity")
  expect_true(all(res_identity$inverse_link_derivatives, na.rm = TRUE), info = "Inverse link derivatives for identity")

  # Mathematical validations for inverse_link
  expect_true(res_inverse$invertibility_theta, info = "Invertibility (theta) for inverse")
  expect_true(res_inverse$invertibility_eta, info = "Invertibility (eta) for inverse")
  expect_true(res_inverse$monotonicity, info = "Monotonicity for inverse")
  expect_true(res_inverse$inverse_theorem, info = "Inverse theorem for inverse")
  expect_true(all(res_inverse$link_derivatives, na.rm = TRUE), info = "Link derivatives for inverse")
  expect_true(all(res_inverse$inverse_link_derivatives, na.rm = TRUE), info = "Inverse link derivatives for inverse")

  # Mathematical validations for inverse_sq_link
  expect_true(res_inverse_sq$invertibility_theta, info = "Invertibility (theta) for inverse_sq")
  # eta domain is strictly positive. Negative eta values result in NaN when passed to linkinv.
  expect_false(res_inverse_sq$invertibility_eta, info = "Invertibility (eta) for inverse_sq (restricted domain)")
  expect_true(res_inverse_sq$monotonicity, info = "Monotonicity for inverse_sq")
  expect_true(res_inverse_sq$inverse_theorem, info = "Inverse theorem for inverse_sq")
  expect_true(all(res_inverse_sq$link_derivatives, na.rm = TRUE), info = "Link derivatives for inverse_sq")
  expect_true(all(res_inverse_sq$inverse_link_derivatives, na.rm = TRUE), info = "Inverse link derivatives for inverse_sq")

  # Mathematical validations for log_link
  expect_true(res_log$invertibility_theta, info = "Invertibility (theta) for log")
  expect_true(res_log$invertibility_eta, info = "Invertibility (eta) for log")
  expect_true(res_log$monotonicity, info = "Monotonicity for log")
  expect_true(res_log$inverse_theorem, info = "Inverse theorem for log")
  expect_true(all(res_log$link_derivatives, na.rm = TRUE), info = "Link derivatives for log")
  expect_true(all(res_log$inverse_link_derivatives, na.rm = TRUE), info = "Inverse link derivatives for log")

  # Mathematical validations for logit_link
  expect_true(res_logit$invertibility_theta, info = "Invertibility (theta) for logit")
  expect_true(res_logit$invertibility_eta, info = "Invertibility (eta) for logit")
  expect_true(res_logit$monotonicity, info = "Monotonicity for logit")
  expect_true(res_logit$inverse_theorem, info = "Inverse theorem for logit")
  expect_true(all(res_logit$link_derivatives, na.rm = TRUE), info = "Link derivatives for logit")
  expect_true(all(res_logit$inverse_link_derivatives, na.rm = TRUE), info = "Inverse link derivatives for logit")

  # Mathematical validations for loglog_link
  expect_true(res_loglog$invertibility_theta, info = "Invertibility (theta) for loglog")
  expect_true(res_loglog$invertibility_eta, info = "Invertibility (eta) for loglog")
  expect_true(res_loglog$monotonicity, info = "Monotonicity for loglog")
  expect_true(res_loglog$inverse_theorem, info = "Inverse theorem for loglog")
  expect_true(all(res_loglog$link_derivatives, na.rm = TRUE), info = "Link derivatives for loglog")
  expect_true(all(res_loglog$inverse_link_derivatives, na.rm = TRUE), info = "Inverse link derivatives for loglog")

  # Mathematical validations for power_link
  expect_true(res_power$invertibility_theta, info = "Invertibility (theta) for power")
  expect_true(res_power$invertibility_eta, info = "Invertibility (eta) for power")
  expect_true(res_power$monotonicity, info = "Monotonicity for power")
  expect_true(res_power$inverse_theorem, info = "Inverse theorem for power")
  expect_true(all(res_power$link_derivatives, na.rm = TRUE), info = "Link derivatives for power")
  expect_true(all(res_power$inverse_link_derivatives, na.rm = TRUE), info = "Inverse link derivatives for power")

  # Mathematical validations for probit_link
  expect_true(res_probit$invertibility_theta, info = "Invertibility (theta) for probit")
  expect_true(res_probit$invertibility_eta, info = "Invertibility (eta) for probit")
  expect_true(res_probit$monotonicity, info = "Monotonicity for probit")
  expect_true(res_probit$inverse_theorem, info = "Inverse theorem for probit")
  expect_true(all(res_probit$link_derivatives, na.rm = TRUE), info = "Link derivatives for probit")
  expect_true(all(res_probit$inverse_link_derivatives, na.rm = TRUE), info = "Inverse link derivatives for probit")

  # Mathematical validations for rhobit_link
  expect_true(res_rhobit$invertibility_theta, info = "Invertibility (theta) for rhobit")
  expect_true(res_rhobit$invertibility_eta, info = "Invertibility (eta) for rhobit")
  expect_true(res_rhobit$monotonicity, info = "Monotonicity for rhobit")
  expect_true(res_rhobit$inverse_theorem, info = "Inverse theorem for rhobit")
  expect_true(all(res_rhobit$link_derivatives, na.rm = TRUE), info = "Link derivatives for rhobit")
  expect_true(all(res_rhobit$inverse_link_derivatives, na.rm = TRUE), info = "Inverse link derivatives for rhobit")

  # Mathematical validations for softplus_link
  expect_true(res_softplus$invertibility_theta, info = "Invertibility (theta) for softplus")
  expect_true(res_softplus$invertibility_eta, info = "Invertibility (eta) for softplus")
  expect_true(res_softplus$monotonicity, info = "Monotonicity for softplus")
  expect_true(res_softplus$inverse_theorem, info = "Inverse theorem for softplus")
  expect_true(all(res_softplus$link_derivatives, na.rm = TRUE), info = "Link derivatives for softplus")
  expect_true(all(res_softplus$inverse_link_derivatives, na.rm = TRUE), info = "Inverse link derivatives for softplus")

  # Mathematical validations for sqrt_link
  expect_true(res_sqrt$invertibility_theta, info = "Invertibility (theta) for sqrt")
  # eta domain is strictly positive. Negative eta evaluates to positive theta, losing the original sign.
  expect_false(res_sqrt$invertibility_eta, info = "Invertibility (eta) for sqrt (restricted domain)")
  expect_true(res_sqrt$monotonicity, info = "Monotonicity for sqrt")
  expect_true(res_sqrt$inverse_theorem, info = "Inverse theorem for sqrt")
  expect_true(all(res_sqrt$link_derivatives, na.rm = TRUE), info = "Link derivatives for sqrt")
  expect_true(all(res_sqrt$inverse_link_derivatives, na.rm = TRUE), info = "Inverse link derivatives for sqrt")

})

# Create a list of all link objects for comprehensive testing of other methods
all_links <- list(
  bounded_double = bounded_link(lwr = 0, upr = 1),
  bounded_lower = bounded_link(lwr = 0),
  cauchit = cauchit_link(),
  cloglog = cloglog_link(),
  identity = identity_link(),
  inverse = inverse_link(),
  inverse_sq = inverse_sq_link(),
  log = log_link(),
  logit = logit_link(),
  loglog = loglog_link(),
  power = power_link(),
  probit = probit_link(),
  rhobit = rhobit_link(),
  softplus = softplus_link(),
  sqrt = sqrt_link()
)

test_that("print.link outputs correct information for all links", {
  for (name in names(all_links)) {
    link_obj <- all_links[[name]]

    # Check that print executes and contains the base string
    expect_output(
      print(link_obj),
      "S7 Link Object:",
      info = paste("print method failed for", name)
    )
  }
})

test_that("plot.link executes without producing errors for all links", {
  # Suppress actual graphic generation to speed up automated testing
  pdf(file = NULL)
  on.exit(dev.off())

  for (name in names(all_links)) {
    link_obj <- all_links[[name]]

    # Expect no errors when plotting within the defined mathematical boundaries
    expect_error(
      plot(link_obj),
      NA,
      info = paste("plot method failed for link:", name)
    )
  }
})

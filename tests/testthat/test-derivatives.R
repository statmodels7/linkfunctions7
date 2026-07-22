# Every link, every derivative order, against an independent reference; plus the
# identities that must hold whatever the link is.

all_links <- function() {
  list(
    identity   = identity_link(),
    log        = log_link(),
    logit      = logit_link(),
    probit     = probit_link(),
    cloglog    = cloglog_link(),
    loglog     = loglog_link(),
    cauchit    = cauchit_link(),
    rhobit     = rhobit_link(),
    sqrt       = sqrt_link(),
    inverse    = inverse_link(),
    inverse_sq = inverse_sq_link(),
    power_2    = power_link(2),
    power_half = power_link(0.5),
    softplus_1 = softplus_link(1),
    softplus_3 = softplus_link(3),
    lower_b    = bounded_link(lwr = 2),
    upper_b    = bounded_link(upr = 5),
    both_b     = bounded_link(lwr = 2, upr = 5)
  )
}

# points well inside the domain, so that numerical differentiation stays there
theta_grid <- function(lk, m = 9) {
  b <- lk@link_bounds
  if (all(is.finite(b))) {
    seq(b[1], b[2], length.out = m + 2)[-c(1, m + 2)]
  } else if (is.finite(b[1])) {
    b[1] + exp(seq(log(0.3), log(4), length.out = m))
  } else if (is.finite(b[2])) {
    b[2] - exp(seq(log(0.3), log(4), length.out = m))
  } else {
    seq(-2.5, 2.5, length.out = m)
  }
}

test_that("analytical derivatives match Richardson extrapolation to fourth order", {
  # The reference differentiates the analytical derivative one order below, so an
  # error at order k cannot hide behind an error at order k-1.
  for (nm in names(all_links())) {
    lk <- all_links()[[nm]]
    th <- theta_grid(lk)
    eta <- linkfun(lk, th)

    for (o in 1:4) {
      ref_f <- vapply(th, function(v) {
        numDeriv::grad(function(z) linkderiv(lk, z, order = o - 1), v)
      }, numeric(1))
      expect_equal(linkderiv(lk, th, order = o), ref_f, tolerance = 1e-5,
        label = sprintf("%s linkderiv order %d", nm, o))

      ref_i <- vapply(eta, function(v) {
        numDeriv::grad(function(z) linkinvderiv(lk, z, order = o - 1), v)
      }, numeric(1))
      expect_equal(linkinvderiv(lk, eta, order = o), ref_i, tolerance = 1e-5,
        label = sprintf("%s linkinvderiv order %d", nm, o))
    }
  }
})

test_that("every link round-trips and satisfies the inverse function theorem", {
  for (nm in names(all_links())) {
    lk <- all_links()[[nm]]
    th <- theta_grid(lk)
    eta <- linkfun(lk, th)

    expect_equal(linkinv(lk, eta), th, tolerance = 1e-10, label = paste(nm, "round-trip"))
    expect_equal(dlinkfun(lk, th) * dlinkinv(lk, eta), rep(1, length(th)),
      tolerance = 1e-10, label = paste(nm, "inverse function theorem"))

    d1 <- dlinkfun(lk, theta_grid(lk, 40))
    expect_true(all(d1 > 0) || all(d1 < 0), label = paste(nm, "strict monotonicity"))
  }
})

test_that("check_link passes on every link the package ships", {
  # It used to report failures for links that are perfectly correct: the
  # evaluation grid came within 1e-3 of the domain boundary, closer than
  # numDeriv's stencil reaches, and the eta-invertibility test was run over a
  # fixed [-4, 4] that a restricted-range link can never produce.
  for (nm in names(all_links())) {
    res <- suppressWarnings(
      utils::capture.output(r <- check_link(all_links()[[nm]]))
    )
    failed <- names(Filter(function(v) !all(v, na.rm = TRUE), r))
    expect_length(failed, 0)
  }
})

test_that("derivatives are vectorised and preserve missing values", {
  # A derivative that reduces to a constant must still return NA for an NA input.
  # R makes this easy to miss: NA^0 is 1, so theta^(lambda - 2) quietly turns a
  # missing parameter into a number as soon as lambda is 2.
  for (nm in names(all_links())) {
    lk <- all_links()[[nm]]
    th <- theta_grid(lk)
    for (o in 0:4) {
      expect_length(linkderiv(lk, th, order = o), length(th))
      expect_length(linkderiv(lk, numeric(0), order = o), 0)
      expect_true(is.na(linkderiv(lk, NA_real_, order = o)),
        label = sprintf("%s linkderiv order %d on NA", nm, o))
      expect_true(is.na(linkinvderiv(lk, NA_real_, order = o)),
        label = sprintf("%s linkinvderiv order %d on NA", nm, o))
    }
  }
})

test_that("an unsupported derivative order is rejected", {
  lk <- log_link()
  expect_error(linkderiv(lk, 1, order = 5), "not supported")
  expect_error(linkinvderiv(lk, 1, order = -1), "not supported")
  expect_error(linkderiv(lk, 1, order = NA), "not supported")
})

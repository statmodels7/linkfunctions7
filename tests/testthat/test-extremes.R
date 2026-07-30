# Behaviour far out in the tails, where the ordinary grids of test-derivatives.R
# never go. Three defects lived here undetected precisely because every other
# test evaluates well inside the domain.

test_that("softplus survives large a * theta", {
  # log(expm1(z)) overflows at z = 709, and the derivatives divide by powers of
  # expm1(z), so the higher orders went first: d4linkfun was NaN from theta = 177
  # with a = 1, and from theta = 24 with a = 30, which is an ordinary value.
  for (a in c(0.5, 1, 30)) {
    lk <- softplus_link(a)
    th <- c(1, 25, 100, 700, 1e4, 1e6)

    expect_true(all(is.finite(linkfun(lk, th))),
      label = sprintf("softplus(a=%g) linkfun finite", a))
    for (o in 1:4) {
      expect_true(all(is.finite(linkderiv(lk, th, order = o))),
        label = sprintf("softplus(a=%g) forward derivative order %d finite", a, o))
    }

    # asymptotically the softplus is the identity, so g(theta) -> theta and
    # g'(theta) -> 1
    big <- 1e4
    expect_equal(linkfun(lk, big), big)
    expect_equal(dlinkfun(lk, big), 1)
    expect_equal(d2linkfun(lk, big), 0)
  }
})

test_that("softplus agrees with the expm1 form wherever that form still works", {
  # The implementation was rewritten in u = -expm1(-a theta); this pins it to the
  # algebraically equivalent expression it replaced, over the range where the old
  # one is still finite.
  for (a in c(0.5, 1, 3)) {
    lk <- softplus_link(a)
    th <- c(0.01, 0.1, 0.5, 1, 3, 10, 30)
    val <- expm1(a * th)

    expect_equal(linkfun(lk, th), log(val) / a)
    expect_equal(dlinkfun(lk, th), (val + 1) / val)
    expect_equal(d2linkfun(lk, th), -a * (val + 1) / val^2)
    expect_equal(d3linkfun(lk, th), a^2 * (val + 1) * (val + 2) / val^3)
    expect_equal(d4linkfun(lk, th),
                 -a^3 * (val + 1) * (val^2 + 6 * val + 6) / val^4)
  }
})

test_that("the exponential links stay exact far into the lower tail", {
  # The floor used to be .Machine$double.eps, which corrupted theta for every
  # eta below about -36 and broke the round trip there.
  eta <- c(-20, -36, -40, -60, -100, -150)

  expect_equal(linkinv(log_link(), eta), exp(eta))
  expect_equal(linkfun(log_link(), linkinv(log_link(), eta)), eta)

  expect_equal(linkinv(bounded_link(lwr = 2), eta) - 2, exp(eta))
  expect_equal(linkinv(cloglog_link(), eta), -expm1(-exp(eta)))
})

test_that("the exponential floor keeps the forward derivatives finite", {
  # The floor is chosen so that -6/theta^4, the binding one, cannot overflow.
  lk <- log_link()
  theta <- linkinv(lk, c(-500, -1000, -Inf))
  for (o in 1:4) {
    expect_true(all(is.finite(linkderiv(lk, theta, order = o))),
      label = sprintf("log link forward derivative order %d finite at the floor", o))
  }
})

test_that("check_link reports a numerical order as numerical, not as passed", {
  # The orders a link does not implement are supplied by the fallbacks, so they
  # can be evaluated -- but checking them would compare a finite difference
  # against a finite difference of the order below, which is the same arithmetic
  # and would agree however wrong the link is. Such an order is therefore left
  # NA, meaning "not checked", and the console says so. Reporting it as PASSED
  # would be the worst of the three possibilities.
  Partial <- S7::new_class("Partial", parent = link)
  S7::method(linkfun,  Partial) <- function(x, theta) log(theta)
  S7::method(linkinv,  Partial) <- function(x, eta)   exp(eta)
  S7::method(dlinkfun, Partial) <- function(x, theta) 1 / theta
  S7::method(dlinkinv, Partial) <- function(x, eta)   exp(eta)

  lk <- Partial(link_name = "partial", link_bounds = c(0, Inf), link_params = NULL)
  txt <- capture.output(out <- suppressWarnings(check_link(lk)))

  expect_equal(attr(out, "analytic_orders"), list(forward = 1L, inverse = 1L))
  expect_true(out$link_derivatives[["order_1"]])
  expect_true(all(is.na(out$link_derivatives[2:4])))
  expect_false(any(out$link_derivatives %in% FALSE))
  expect_match(paste(txt, collapse = " "), "numerical")

  # a link implementing nothing at all is reported as wholly numerical
  Bare <- S7::new_class("Bare", parent = link)
  S7::method(linkfun, Bare) <- function(x, theta) log(theta)
  S7::method(linkinv, Bare) <- function(x, eta)   exp(eta)
  bare <- Bare(link_name = "bare", link_bounds = c(0, Inf), link_params = NULL)
  txt2 <- capture.output(out2 <- suppressWarnings(check_link(bare)))
  expect_equal(attr(out2, "analytic_orders"), list(forward = 0L, inverse = 0L))
  expect_true(all(is.na(out2$link_derivatives)))
})

test_that("the numerical fallbacks reproduce a known link", {
  # The log link, defined without any of its derivatives, against the closed
  # forms. Accuracy degrades with order, as a stencil must, but every order is
  # usable -- and supplying the analytical second derivative makes the third and
  # fourth far better, because the fallback then differentiates once instead of
  # standing on three differences of its own.
  Bare <- S7::new_class("BareLog", parent = link)
  S7::method(linkfun, Bare) <- function(x, theta) log(theta)
  S7::method(linkinv, Bare) <- function(x, eta)   exp(eta)
  bare <- Bare(link_name = "bare log", link_bounds = c(0, Inf), link_params = NULL)

  Half <- S7::new_class("HalfLog", parent = link)
  S7::method(linkfun,   Half) <- function(x, theta) log(theta)
  S7::method(linkinv,   Half) <- function(x, eta)   exp(eta)
  S7::method(dlinkfun,  Half) <- function(x, theta) 1 / theta
  S7::method(d2linkfun, Half) <- function(x, theta) -1 / theta^2
  S7::method(dlinkinv,  Half) <- function(x, eta) exp(eta)
  S7::method(d2linkinv, Half) <- function(x, eta) exp(eta)
  half <- Half(link_name = "half log", link_bounds = c(0, Inf), link_params = NULL)

  th <- c(0.5, 1, 2, 5)
  exact <- list(1 / th, -1 / th^2, 2 / th^3, -6 / th^4)
  tol <- c(1e-8, 1e-6, 1e-4, 1e-3)          # what a stencil of that order can do

  for (k in 1:4) {
    e_bare <- max(abs(linkderiv(bare, th, order = k) - exact[[k]]) / abs(exact[[k]]))
    e_half <- max(abs(linkderiv(half, th, order = k) - exact[[k]]) / abs(exact[[k]]))
    expect_lt(e_bare, tol[k])
    expect_lt(e_half, tol[k])
  }

  # the payoff of not nesting: at fourth order, knowing two orders analytically
  # is worth orders of magnitude
  e4_bare <- max(abs(linkderiv(bare, th, order = 4) - exact[[4]]) / abs(exact[[4]]))
  e4_half <- max(abs(linkderiv(half, th, order = 4) - exact[[4]]) / abs(exact[[4]]))
  expect_lt(e4_half, e4_bare / 100)

  # and the inverse direction, where every derivative of exp is exp
  expect_equal(linkinvderiv(bare, 0.7, order = 3), exp(0.7), tolerance = 1e-4)
})

test_that("the shipped links are all analytic to fourth order", {
  # The fallbacks exist for user-defined links; nothing in the catalogue should
  # ever reach them.
  for (nm in names(all_links())) {
    expect_equal(link_fallback_orders(all_links()[[nm]]),
                 list(forward = 4L, inverse = 4L),
                 label = nm)
  }
})

test_that("a fallback respects the domain and propagates missingness", {
  Bare <- S7::new_class("BareLog2", parent = link)
  S7::method(linkfun, Bare) <- function(x, theta) log(theta)
  S7::method(linkinv, Bare) <- function(x, eta)   exp(eta)
  bare <- Bare(link_name = "bare log", link_bounds = c(0, Inf), link_params = NULL)

  # close to the boundary the step must shrink rather than step through zero
  expect_true(all(is.finite(linkderiv(bare, 1e-6, order = 4))))
  expect_true(is.na(linkderiv(bare, NA_real_, order = 2)))
  expect_length(linkderiv(bare, numeric(0), order = 3), 0)
})

test_that("check_link still catches a derivative that is merely wrong", {
  # The guard above must not have been bought by weakening the check.
  Wrong <- S7::new_class("Wrong", parent = link)
  S7::method(linkfun,   Wrong) <- function(x, theta) log(theta)
  S7::method(linkinv,   Wrong) <- function(x, eta)   exp(eta)
  S7::method(dlinkfun,  Wrong) <- function(x, theta) 1 / theta
  S7::method(d2linkfun, Wrong) <- function(x, theta) -1.05 / theta^2  # 5% out
  S7::method(d3linkfun, Wrong) <- function(x, theta) 2 / theta^3
  S7::method(d4linkfun, Wrong) <- function(x, theta) -6 / theta^4
  S7::method(dlinkinv,  Wrong) <- function(x, eta) exp(eta)
  S7::method(d2linkinv, Wrong) <- function(x, eta) exp(eta)
  S7::method(d3linkinv, Wrong) <- function(x, eta) exp(eta)
  S7::method(d4linkinv, Wrong) <- function(x, eta) exp(eta)

  lk <- Wrong(link_name = "wrong", link_bounds = c(0, Inf), link_params = NULL)
  res <- suppressWarnings(capture.output(out <- check_link(lk)))

  expect_true(out$link_derivatives[["order_1"]])
  expect_false(out$link_derivatives[["order_2"]])
  # and an unbroken link of the same shape must still pass, so the check cannot
  # be failing trivially
  res2 <- suppressWarnings(capture.output(ok <- check_link(log_link())))
  expect_true(all(unlist(ok)))
})

test_that("the shared logistic polynomials match their expanded forms", {
  # logit, doubly bounded and softplus all route through logistic_deriv(); the
  # Horner evaluation must agree with the expanded polynomials it replaced.
  eta <- seq(-8, 8, length.out = 101)
  p <- stats::plogis(eta)

  expect_equal(dlinkinv(logit_link(), eta), p * (1 - p))
  expect_equal(d2linkinv(logit_link(), eta), p * (1 - p) * (1 - 2 * p))
  expect_equal(d3linkinv(logit_link(), eta), p * (1 - p) * (1 - 6 * p + 6 * p^2))
  expect_equal(d4linkinv(logit_link(), eta),
               p * (1 - p) * (1 - 14 * p + 36 * p^2 - 24 * p^3))

  db <- bounded_link(-3, 7)
  expect_equal(db@width, 10)
  expect_equal(d4linkinv(db, eta),
               10 * p * (1 - p) * (1 - 14 * p + 36 * p^2 - 24 * p^3))

  a <- 2
  sp <- softplus_link(a)
  q <- stats::plogis(a * eta)
  expect_equal(d2linkinv(sp, eta), a * q * (1 - q))
  expect_equal(d4linkinv(sp, eta), a^3 * q * (1 - q) * (1 - 6 * q + 6 * q^2))
})

test_that("constructors reject malformed arguments", {
  expect_error(bounded_link(lwr = c(1, 2)), "single finite number")
  expect_error(bounded_link(lwr = NA), "single finite number")
  expect_error(bounded_link(upr = Inf), "single finite number")
  expect_error(softplus_link(a = 0), "strictly greater than 0")
  expect_error(softplus_link(a = -1), "strictly greater than 0")
  expect_error(softplus_link(a = c(1, 2)), "single number")
  expect_error(power_link(lambda = NA), "single finite number")
  expect_error(power_link(lambda = c(1, 2)), "single finite number")
})

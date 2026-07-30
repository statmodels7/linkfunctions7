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

test_that("check_link fails a link whose derivatives are not implemented", {
  # An unimplemented derivative used to leave NA, and the summary reduced with
  # na.rm = TRUE, so a link supplying only its first derivative reported PASSED
  # on all four orders.
  Partial <- S7::new_class("Partial", parent = link)
  S7::method(linkfun,  Partial) <- function(x, theta) log(theta)
  S7::method(linkinv,  Partial) <- function(x, eta)   exp(eta)
  S7::method(dlinkfun, Partial) <- function(x, theta) 1 / theta
  S7::method(dlinkinv, Partial) <- function(x, eta)   exp(eta)

  lk <- Partial(link_name = "partial", link_bounds = c(0, Inf), link_params = NULL)
  res <- suppressWarnings(capture.output(out <- check_link(lk)))

  expect_true(out$link_derivatives[["order_1"]])
  expect_false(any(out$link_derivatives[2:4]))
  expect_false(anyNA(out$link_derivatives))
  expect_false(anyNA(out$inverse_link_derivatives))
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

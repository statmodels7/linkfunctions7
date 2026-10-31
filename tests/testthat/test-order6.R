# The sixth derivative of every link, in both directions.
#
# Stage one of the order: the generics, the routers, the numerical fallback of
# the base class and the `case 6:` slots of the compiled kernels. A link that
# registers no sixth-order method is answered by one central first-difference
# stencil applied to its analytic fifth derivative, never a chain.

order6_grid <- function(lk, m = 9) {
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

test_that("every shipped link answers order six in both directions", {
  for (nm in names(all_links())) {
    lk <- all_links()[[nm]]
    th <- order6_grid(lk)
    eta <- linkfun(lk, th)
    f <- d6linkfun(lk, th)
    i <- d6linkinv(lk, eta)
    expect_length(f, length(th))
    expect_length(i, length(eta))
    expect_true(!anyNA(f), label = paste(nm, "d6linkfun"))
    expect_true(!anyNA(i), label = paste(nm, "d6linkinv"))
  }
})

test_that("the order-six answer matches Richardson on the analytic fifth", {
  # one Richardson pass on the analytic order 5, the rule of check_link();
  # compared on max(1, |reference|) because several sixth derivatives cross
  # zero inside the grid
  for (nm in names(all_links())) {
    lk <- all_links()[[nm]]
    th <- order6_grid(lk)
    eta <- linkfun(lk, th)
    ref_f <- vapply(th, function(v)
      numDeriv::grad(function(z) d5linkfun(lk, z), v), numeric(1))
    ref_i <- vapply(eta, function(v)
      numDeriv::grad(function(z) d5linkinv(lk, z), v), numeric(1))
    ef <- max(abs(d6linkfun(lk, th) - ref_f) / pmax(1, abs(ref_f)))
    ei <- max(abs(d6linkinv(lk, eta) - ref_i) / pmax(1, abs(ref_i)))
    expect_lt(ef, 1e-5, label = paste(nm, "d6linkfun"))
    expect_lt(ei, 1e-5, label = paste(nm, "d6linkinv"))
  }
})

test_that("the fallback is one stencil on the analytic fifth", {
  # the log link, whose sixth derivatives are known in closed form
  Bare <- S7::new_class("BareLog6", parent = link)
  S7::method(linkfun, Bare) <- function(x, theta) log(theta)
  S7::method(linkinv, Bare) <- function(x, eta) exp(eta)
  bare <- Bare(link_name = "bare log", link_bounds = c(0, Inf), link_params = NULL)
  expect_equal(link_fallback_orders(bare), list(forward = 0L, inverse = 0L))
  expect_equal(d6linkinv(bare, 0.3), exp(0.3), tolerance = 1e-2)

  Half <- S7::new_class("HalfLog6", parent = link)
  S7::method(linkfun, Half) <- function(x, theta) log(theta)
  S7::method(linkinv, Half) <- function(x, eta) exp(eta)
  S7::method(dlinkfun, Half) <- function(x, theta) 1 / theta
  S7::method(d2linkfun, Half) <- function(x, theta) -1 / theta^2
  S7::method(d3linkfun, Half) <- function(x, theta) 2 / theta^3
  S7::method(d4linkfun, Half) <- function(x, theta) -6 / theta^4
  S7::method(d5linkfun, Half) <- function(x, theta) 24 / theta^5
  S7::method(dlinkinv, Half) <- function(x, eta) exp(eta)
  S7::method(d2linkinv, Half) <- function(x, eta) exp(eta)
  S7::method(d3linkinv, Half) <- function(x, eta) exp(eta)
  S7::method(d4linkinv, Half) <- function(x, eta) exp(eta)
  S7::method(d5linkinv, Half) <- function(x, eta) exp(eta)
  half <- Half(link_name = "half log", link_bounds = c(0, Inf), link_params = NULL)
  expect_equal(link_fallback_orders(half), list(forward = 5L, inverse = 5L))
  th <- c(0.5, 1, 2, 7)
  expect_equal(d6linkfun(half, th), -120 / th^6, tolerance = 1e-6)
  expect_equal(d6linkinv(half, c(-1, 0, 2)), exp(c(-1, 0, 2)), tolerance = 1e-7)

  # check_link leaves the order the fallback supplies unchecked
  out <- NULL
  txt <- capture.output(out <- check_link(half))
  expect_true(is.na(out$link_derivatives[["order_6"]]))
  expect_true(is.na(out$inverse_link_derivatives[["order_6"]]))
  expect_true(all(out$link_derivatives[1:5]))
  expect_true(any(grepl("PASSED to order 5, 1 numerical", txt)))
})

test_that("the routers reach order six and refuse the seventh", {
  expect_equal(linkderiv(log_link(), 2, order = 6), d6linkfun(log_link(), 2))
  expect_equal(linkinvderiv(log_link(), 1, order = 6), d6linkinv(log_link(), 1))
  expect_error(linkderiv(log_link(), 2, order = 7), "not supported")
  expect_error(linkinvderiv(log_link(), 2, order = 7), "not supported")
  # linkinv_derivs() carries the sixth
  l6 <- linkinv_derivs(logit_link(), 0.4, order = 6)
  expect_length(l6, 6L)
  expect_equal(l6[[6]], d6linkinv(logit_link(), 0.4))
  expect_length(linkinv_derivs(log_link(), 0.4, order = 6), 6L)
})

test_that("a compiled kernel answers NA at order six until it carries it", {
  # pins the `case 6:` slots: an unsupported order is missing, never the
  # value of another order. Stage two replaces these with closed forms, and
  # this test then moves to order 7.
  e <- c(-1.2, 0.3, 1.4)
  p <- stats::plogis(e)
  th <- c(0.2, 0.5, 0.8)
  expect_true(all(is.na(lk_logit_fwd_cpp(th, 7L))))
  expect_true(all(is.na(lk_probit_fwd_cpp(th, 7L))))
  expect_true(all(is.na(lk_cloglog_fwd_cpp(th, 7L))))
  expect_true(all(is.na(lk_loglog_fwd_cpp(th, 7L))))
  expect_true(all(is.na(lk_cauchit_fwd_cpp(th, 7L))))
  expect_true(all(is.na(lk_rhobit_fwd_cpp(th, 7L))))
  expect_true(all(is.na(lk_softplus_fwd_cpp(th, 1, 7L))))
  expect_true(all(is.na(lk_probit_inv_cpp(e, 7L))))
  expect_true(all(is.na(lk_cloglog_inv_cpp(e, 7L))))
  expect_true(all(is.na(lk_loglog_inv_cpp(e, 7L))))
  expect_true(all(is.na(lk_cauchit_inv_cpp(e, 7L))))
  expect_true(all(is.na(lk_rhobit_inv_cpp(e, 7L))))
  expect_true(all(is.na(lk_logit_inv_cpp(e, 7L))))
  expect_true(all(is.na(lk_logistic_poly_cpp(p, 7L))))
  # stage-one slots, to be replaced by closed forms
  expect_true(all(is.na(lk_logit_inv_cpp(e, 6L))))
  expect_true(all(is.na(lk_probit_fwd_cpp(th, 6L))))
})

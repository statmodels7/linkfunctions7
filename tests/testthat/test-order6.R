# The sixth derivative of every link, in both directions.
#
# The sixth derivative of every link. A link that registers no sixth-order
# method is answered by one central first-difference stencil applied to its
# analytic fifth derivative, never a chain; the shipped links are all analytic.

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

test_that("a compiled kernel answers NA above order six", {
  # an unsupported order is missing, never the value of another order
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
  # and order six is carried
  expect_true(all(is.finite(lk_logit_inv_cpp(e, 6L))))
  expect_true(all(is.finite(lk_probit_fwd_cpp(th, 6L))))
})

test_that("every shipped link is analytic to order six, and check_link says so", {
  for (nm in names(all_links())) {
    lk <- all_links()[[nm]]
    expect_equal(link_fallback_orders(lk), list(forward = 6L, inverse = 6L), label = nm)
    txt <- capture.output(out <- check_link(lk))
    expect_true(all(unlist(out)), label = paste(nm, "check_link"))
    expect_equal(sum(grepl("[PASSED]", txt, fixed = TRUE)), 6L, label = nm)
    expect_false(any(grepl("numerical", txt)), label = nm)
  }
})

test_that("order six matches Richardson on the analytic fifth over a wide grid", {
  # a wider grid than above, with the tails a model reaches, one Richardson
  # pass on the analytic fifth; relative to max(1, |reference|)
  wide <- function(lk) {
    b <- lk@link_bounds
    if (all(is.finite(b))) {
      u <- c(1e-3, 1e-2, seq(0.05, 0.95, length.out = 19), 0.99, 1 - 1e-3)
      b[1] + u * (b[2] - b[1])
    } else if (is.finite(b[1])) {
      b[1] + exp(seq(log(0.05), log(30), length.out = 25))
    } else if (is.finite(b[2])) {
      b[2] - exp(seq(log(0.05), log(30), length.out = 25))
    } else seq(-6, 6, length.out = 25)
  }
  for (nm in names(all_links())) {
    lk <- all_links()[[nm]]
    th <- wide(lk)
    eta <- linkfun(lk, th)
    ref_f <- vapply(th, function(v)
      numDeriv::grad(function(z) d5linkfun(lk, z), v), numeric(1))
    ref_i <- vapply(eta, function(v)
      numDeriv::grad(function(z) d5linkinv(lk, z), v), numeric(1))
    ef <- max(abs(d6linkfun(lk, th) - ref_f) / pmax(1, abs(ref_f)))
    ei <- max(abs(d6linkinv(lk, eta) - ref_i) / pmax(1, abs(ref_i)))
    expect_lt(ef, 1e-6, label = paste(nm, "d6linkfun"))
    expect_lt(ei, 1e-6, label = paste(nm, "d6linkinv"))
  }
})

test_that("a sixth derivative with a coefficient 5 per cent out is caught", {
  # negative control of the test above: the logit's sixth inverse derivative
  # with its leading polynomial coefficient perturbed must read a large error
  eta <- seq(-4, 4, length.out = 21)
  p <- stats::plogis(eta); pq <- p * (1 - p)
  good <- pq * (1 - 2 * p) * (1 + pq * (-60 + 360 * pq))
  bad <- pq * (1 - 2 * p) * (1 + pq * (-60 * 1.05 + 360 * pq))
  expect_equal(d6linkinv(logit_link(), eta), good, tolerance = 1e-12)
  expect_gt(max(abs(d6linkinv(logit_link(), eta) - bad)), 1e-4)
})

test_that("the sixth derivatives follow the inverse function theorem", {
  # g(h(eta)) = eta differentiated six times: the sixth inverse derivative is
  # fixed by the forward derivatives and the inverse's first five, a route
  # that shares no arithmetic with the kernels. Checked as a residual of the
  # relation sum_j g^(j)(h) B_{6,j}(h', ..., h^(6-j+1)) = 0, using the partial
  # Bell polynomials with the sixth inverse derivative left free.
  for (nm in setdiff(names(all_links()), c("identity", "sqrt", "power_2"))) {
    lk <- all_links()[[nm]]
    b <- lk@link_bounds
    th <- if (all(is.finite(b))) b[1] + c(0.2, 0.4, 0.6, 0.8) * (b[2] - b[1]) else
      if (is.finite(b[1])) b[1] + c(0.5, 1, 2, 4) else
      if (is.finite(b[2])) b[2] - c(0.5, 1, 2, 4) else c(-1.5, -0.3, 0.4, 1.7)
    eta <- linkfun(lk, th)
    h <- lapply(list(dlinkinv, d2linkinv, d3linkinv, d4linkinv, d5linkinv, d6linkinv),
                function(f) f(lk, eta))
    g <- lapply(list(dlinkfun, d2linkfun, d3linkfun, d4linkfun, d5linkfun, d6linkfun),
                function(f) f(lk, th))
    terms <- list(
      g[[1]] * h[[6]],
      g[[2]] * (6 * h[[1]] * h[[5]] + 15 * h[[2]] * h[[4]] + 10 * h[[3]]^2),
      g[[3]] * (15 * h[[1]]^2 * h[[4]] + 60 * h[[1]] * h[[2]] * h[[3]] + 15 * h[[2]]^3),
      g[[4]] * (20 * h[[1]]^3 * h[[3]] + 45 * h[[1]]^2 * h[[2]]^2),
      g[[5]] * 15 * h[[1]]^4 * h[[2]],
      g[[6]] * h[[1]]^6)
    res <- Reduce(`+`, terms)
    scale <- Reduce(`+`, lapply(terms, abs))
    expect_lt(max(abs(res) / scale), 1e-9, label = nm)
  }
})

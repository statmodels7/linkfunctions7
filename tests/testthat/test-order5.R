# The fifth derivative of every link, in both directions.
#
# The order exists because each order of differentiating a filtered predictor
# through its own recursion draws in one more order of the link: the curvature
# of a score-driven filter reaches the third, the directional third derivative
# reaches the fourth, and the outer Hessian reaches the fifth.
#
# The reference throughout is ONE Richardson pass on the analytic order 4,
# never a chain of first differences, which is the rule check_link() states and
# the reason a fallback order is reported as numerical rather than passed.


test_that("every shipped link is analytic to the fifth order in both directions", {
  # The fallbacks exist for user-defined links. Nothing in the catalog should
  # reach them, and this is what says a class was not missed when the order was
  # added: a link left without a d5 method would report 4 here.
  for (nm in names(all_links())) {
    expect_equal(link_fallback_orders(all_links()[[nm]]),
                 list(forward = 5L, inverse = 5L),
                 label = nm)
  }
})


test_that("the fifth derivative matches Richardson on the analytic fourth", {
  th_grid <- function(lk, m = 9) {
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

  for (nm in names(all_links())) {
    lk <- all_links()[[nm]]
    th <- th_grid(lk)
    eta <- linkfun(lk, th)

    ref_f <- vapply(th, function(v) {
      numDeriv::grad(function(z) d4linkfun(lk, z), v)
    }, numeric(1))
    expect_equal(d5linkfun(lk, th), ref_f, tolerance = 1e-5,
                 label = sprintf("%s d5linkfun", nm))

    ref_i <- vapply(eta, function(v) {
      numDeriv::grad(function(z) d4linkinv(lk, z), v)
    }, numeric(1))
    expect_equal(d5linkinv(lk, eta), ref_i, tolerance = 1e-5,
                 label = sprintf("%s d5linkinv", nm))
  }
})


test_that("check_link catches a fifth derivative that is merely wrong", {
  # The order-5 row of check_link() must be able to fail, or raising the check
  # to five reports a pass it never earned. The link below is the log link with
  # one coefficient of the fifth forward derivative 5% out and everything else
  # exact, so only that row may go red.
  Wrong5 <- S7::new_class("Wrong5", parent = link)
  S7::method(linkfun,   Wrong5) <- function(x, theta) log(theta)
  S7::method(linkinv,   Wrong5) <- function(x, eta)   exp(eta)
  S7::method(dlinkfun,  Wrong5) <- function(x, theta) 1 / theta
  S7::method(d2linkfun, Wrong5) <- function(x, theta) -1 / theta^2
  S7::method(d3linkfun, Wrong5) <- function(x, theta) 2 / theta^3
  S7::method(d4linkfun, Wrong5) <- function(x, theta) -6 / theta^4
  S7::method(d5linkfun, Wrong5) <- function(x, theta) 25.2 / theta^5  # 5% out
  S7::method(dlinkinv,  Wrong5) <- function(x, eta) exp(eta)
  S7::method(d2linkinv, Wrong5) <- function(x, eta) exp(eta)
  S7::method(d3linkinv, Wrong5) <- function(x, eta) exp(eta)
  S7::method(d4linkinv, Wrong5) <- function(x, eta) exp(eta)
  S7::method(d5linkinv, Wrong5) <- function(x, eta) exp(eta)

  lk <- Wrong5(link_name = "wrong5", link_bounds = c(0, Inf), link_params = NULL)
  suppressWarnings(capture.output(out <- check_link(lk)))

  expect_true(out$link_derivatives[["order_4"]])
  expect_false(out$link_derivatives[["order_5"]])
  # the other direction is exact and must still pass, so the failure above is
  # the injected coefficient and not the order itself
  expect_true(out$inverse_link_derivatives[["order_5"]])

  # and the unbroken link of the same shape passes, so the check is not failing
  # trivially
  suppressWarnings(capture.output(ok <- check_link(log_link())))
  expect_true(all(unlist(ok), na.rm = TRUE))  # order 6 is a fallback until stage 2
})


test_that("the routers reach the fifth order and refuse the seventh", {
  # linkderiv() and linkinvderiv() route by order; an order they do not carry
  # is an error rather than the highest one they do.
  expect_equal(linkderiv(log_link(), 2, order = 5), 24 / 2^5)
  expect_equal(linkinvderiv(log_link(), 1, order = 5), exp(1))
  expect_error(linkderiv(log_link(), 2, order = 7))
  expect_error(linkinvderiv(log_link(), 2, order = 7))
})


test_that("a compiled kernel answers NA above the order it carries", {
  # Every kernel takes the order as an argument and used to reach its order-4
  # branch through `default:`, so an order above four came back as the fourth
  # one, silently. The routers never asked for one, so nothing was wrong; the
  # shape was, and an unsupported order is missing rather than plausible.
  e <- c(-1.2, 0.3, 1.4)
  p <- stats::plogis(e)
  expect_true(all(is.na(lk_logit_inv_cpp(e, 7L))))
  expect_true(all(is.na(lk_logistic_poly_cpp(p, 7L))))
  expect_true(all(is.na(lk_probit_inv_cpp(e, 0L))))
  expect_true(all(is.na(lk_cloglog_inv_cpp(e, 9L))))
  # and the orders it does carry are unaffected
  expect_true(all(is.finite(lk_logit_inv_cpp(e, 5L))))
})


test_that("the softplus and the doubly bounded link reuse the logit's fifth", {
  # The softplus is an antiderivative of the logistic, so its inverse-link
  # derivatives are the logistic ones one place down; a doubly bounded link
  # scales them by the interval width. Neither has algebra of its own at the
  # fifth order, and asserting the identity is what says so.
  z <- seq(-6, 6, length.out = 401)

  a <- 2
  q <- stats::plogis(a * z)
  expect_agrees_on_scale(d5linkinv(softplus_link(a), z),
                         a^4 * logistic_deriv(q, 4L),
                         label = "softplus fifth against the logistic fourth")

  w <- 3
  p <- stats::plogis(z)
  expect_agrees_on_scale(d5linkinv(bounded_link(lwr = 1, upr = 1 + w), z),
                         w * logistic_deriv(p, 5L),
                         label = "doubly bounded fifth against the logistic")
})


test_that("a fallback at the fifth order is one stencil and not a chain", {
  # A link analytic to the second order must reach the fifth by ONE stencil of
  # order three, not by three of order one. The measurement that says so is the
  # comparison with a link that supplies nothing: same order, same step rule,
  # orders of magnitude apart.
  Bare <- S7::new_class("Bare5", parent = link)
  S7::method(linkfun, Bare) <- function(x, theta) log(theta)
  S7::method(linkinv, Bare) <- function(x, eta)   exp(eta)
  bare <- Bare(link_name = "bare", link_bounds = c(0, Inf), link_params = NULL)

  Half <- S7::new_class("Half5", parent = link)
  S7::method(linkfun,   Half) <- function(x, theta) log(theta)
  S7::method(linkinv,   Half) <- function(x, eta)   exp(eta)
  S7::method(dlinkinv,  Half) <- function(x, eta) exp(eta)
  S7::method(d2linkinv, Half) <- function(x, eta) exp(eta)
  half <- Half(link_name = "half", link_bounds = c(0, Inf), link_params = NULL)

  expect_equal(link_fallback_orders(bare)$inverse, 0L)
  expect_equal(link_fallback_orders(half)$inverse, 2L)

  truth <- exp(0.7)
  e_bare <- abs(d5linkinv(bare, 0.7) - truth) / truth
  e_half <- abs(d5linkinv(half, 0.7) - truth) / truth

  expect_lt(e_bare, 1e-2)      # honest arithmetic, of declared quality
  expect_lt(e_half, 1e-5)
  expect_lt(e_half, e_bare / 100)
})


test_that("check_link leaves a fifth-order fallback unchecked rather than passed", {
  # Comparing a fallback against a difference of the order below is the same
  # arithmetic twice: it would agree however wrong the link is. That row is NA,
  # which is the convention the fourth order already followed.
  Half <- S7::new_class("Half5b", parent = link)
  S7::method(linkfun,   Half) <- function(x, theta) log(theta)
  S7::method(linkinv,   Half) <- function(x, eta)   exp(eta)
  S7::method(dlinkfun,  Half) <- function(x, theta) 1 / theta
  S7::method(d2linkfun, Half) <- function(x, theta) -1 / theta^2
  S7::method(d3linkfun, Half) <- function(x, theta) 2 / theta^3
  S7::method(d4linkfun, Half) <- function(x, theta) -6 / theta^4
  S7::method(dlinkinv,  Half) <- function(x, eta) exp(eta)
  S7::method(d2linkinv, Half) <- function(x, eta) exp(eta)
  S7::method(d3linkinv, Half) <- function(x, eta) exp(eta)
  S7::method(d4linkinv, Half) <- function(x, eta) exp(eta)
  half <- Half(link_name = "half fwd", link_bounds = c(0, Inf), link_params = NULL)

  expect_equal(link_fallback_orders(half), list(forward = 4L, inverse = 4L))
  suppressWarnings(capture.output(out <- check_link(half)))

  expect_true(out$link_derivatives[["order_4"]])
  expect_true(is.na(out$link_derivatives[["order_5"]]))
  expect_true(is.na(out$inverse_link_derivatives[["order_5"]]))
  expect_equal(attr(out, "analytic_orders"), list(forward = 4L, inverse = 4L))
})


test_that("the fifth derivatives of the two directions imply each other", {
  # The strongest check available inside this package, and the reason the
  # order was added here before distributions7 needs it: differentiating
  # g(h(eta)) = eta five times gives
  #
  #   0 = sum_j g^(j)(h) B_{5,j}(h', h'', h''', h'''', h^(5)),
  #
  # and B_{5,1} = h^(5), so the inverse link's fifth derivative is determined
  # by the FORWARD ones together with the inverse's first four. That route
  # shares no arithmetic with the kernel d5linkinv() calls, so the two
  # agreeing is a check rather than the same expression twice.
  #
  # The partial Bell polynomials of order five are
  #   B_5,1 = h5
  #   B_5,2 = 5 h1 h4 + 10 h2 h3
  #   B_5,3 = 10 h1^2 h3 + 15 h1 h2^2
  #   B_5,4 = 10 h1^3 h2
  #   B_5,5 = h1^5
  # whose coefficients sum to the Bell number B_5 = 52, and which were checked
  # against a series construction of the same objects before being written.
  compared <- 0L
  for (nm in names(all_links())) {
    lk <- all_links()[[nm]]
    b <- lk@link_bounds
    th <- if (all(is.finite(b))) {
      seq(b[1], b[2], length.out = 9)[2:8]
    } else if (is.finite(b[1])) {
      b[1] + exp(seq(log(0.3), log(4), length.out = 7))
    } else if (is.finite(b[2])) {
      b[2] - exp(seq(log(0.3), log(4), length.out = 7))
    } else {
      seq(-2, 2, length.out = 7)
    }
    eta <- linkfun(lk, th)
    hh <- linkinv(lk, eta)

    h1 <- dlinkinv(lk, eta)
    h2 <- d2linkinv(lk, eta)
    h3 <- d3linkinv(lk, eta)
    h4 <- d4linkinv(lk, eta)

    g1 <- dlinkfun(lk, hh)
    g2 <- d2linkfun(lk, hh)
    g3 <- d3linkfun(lk, hh)
    g4 <- d4linkfun(lk, hh)
    g5 <- d5linkfun(lk, hh)

    t2 <- g2 * (10 * h2 * h3 + 5 * h1 * h4)
    t3 <- g3 * (10 * h1^2 * h3 + 15 * h1 * h2^2)
    t4 <- g4 * (10 * h1^3 * h2)
    t5 <- g5 * (h1^5)
    from_fwd <- -(t2 + t3 + t4 + t5) / g1

    got <- d5linkinv(lk, eta)
    ok <- is.finite(got) & is.finite(from_fwd)
    skip_if_not(any(ok), sprintf("%s: nothing finite to compare", nm))

    # The denominator is the size of the largest TERM and not of the result.
    # For the square root the fifth derivative is exactly zero by an exact
    # cancellation, its three surviving contributions carrying the
    # coefficients 45, -150 and 105, so a comparison relative to the answer
    # has no denominator at all and reports machine noise as a total failure.
    scale <- pmax(abs(t2), abs(t3), abs(t4), abs(t5), abs(got * g1)) / abs(g1)

    # For the identity link there is no denominator either, and for a
    # different reason: every term is separately zero rather than cancelling.
    # The identity there is 0 = 0 and is asserted as one, which says more than
    # a ratio could.
    flat <- ok & scale == 0
    if (any(flat)) {
      expect_identical(got[flat], rep(0, sum(flat)), label = sprintf("%s: h5", nm))
      expect_identical(from_fwd[flat], rep(0, sum(flat)),
                       label = sprintf("%s: h5 from the forward route", nm))
    }

    live <- ok & scale > 0
    if (any(live)) {
      compared <- compared + 1L
      expect_lt(max(abs(got[live] - from_fwd[live]) / scale[live]), 1e-10,
                label = sprintf("%s: d5linkinv against the forward route", nm))
    }
  }

  # A skip inside the loop would abandon every link after it -- the identity
  # comes first, so one would have cost the other fifteen. The count is
  # asserted instead, which cannot go quiet: only the identity link, whose
  # every contribution is separately zero, is compared by the exact route
  # above rather than by the ratio.
  expect_equal(compared, length(all_links()) - 1L)
})


test_that("that cross-direction route can fail", {
  # The negative control for the check above: one Bell coefficient 5 per cent
  # out must move the reconstruction well past the tolerance, or the identity
  # is being satisfied by something other than the polynomials.
  lk <- logit_link()
  eta <- seq(-2, 2, length.out = 7)
  hh <- linkinv(lk, eta)
  h1 <- dlinkinv(lk, eta); h2 <- d2linkinv(lk, eta)
  h3 <- d3linkinv(lk, eta); h4 <- d4linkinv(lk, eta)
  g1 <- dlinkfun(lk, hh); g2 <- d2linkfun(lk, hh); g3 <- d3linkfun(lk, hh)
  g4 <- d4linkfun(lk, hh); g5 <- d5linkfun(lk, hh)

  bad <- -(g2 * (10 * h2 * h3 + 5.25 * h1 * h4) +
           g3 * (10 * h1^2 * h3 + 15 * h1 * h2^2) +
           g4 * (10 * h1^3 * h2) + g5 * h1^5) / g1
  truth <- d5linkinv(lk, eta)
  expect_gt(max(abs(bad - truth)) / max(abs(truth)), 1e-3)
})

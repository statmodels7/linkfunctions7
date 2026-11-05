# The scalar C entry points of the fast route: each one mirrors the R
# method it stands for, and the comparison is identical(), not a tolerance
# -- the consumer's twin test rests on these being the same numbers to the
# bit.

test_that("the scalar entries mirror the R methods bit for bit", {
  eta <- c(-800, -178, -40, -1, 0, 0.5, 3, 37, 700, 710)
  cases <- list(list(l = identity_link(), n = "identity"),
                list(l = log_link(), n = "log"))
  for (cs in cases) {
    pr <- lf7_scalar_probe(cs$n, eta, as.numeric(cs$l@link_bounds))
    expect_gte(pr$id, 0)
    m_inv <- S7::method(linkinv, S7::S7_class(cs$l))
    m_d1 <- S7::method(dlinkinv, S7::S7_class(cs$l))
    m_d2 <- S7::method(d2linkinv, S7::S7_class(cs$l))
    expect_identical(pr$h, as.numeric(m_inv(cs$l, eta)))
    expect_identical(pr$h1, as.numeric(m_d1(cs$l, eta)))
    expect_identical(pr$h2, as.numeric(m_d2(cs$l, eta)))
    expect_identical(pr$clamped,
                     as.numeric(link_bounds_clamp(m_inv(cs$l, eta),
                                                  cs$l@link_bounds)))
  }
})

test_that("an unknown link answers -1 and the probe says so", {
  pr <- lf7_scalar_probe("no-such-link", c(0, 1), c(-Inf, Inf))
  expect_identical(pr$id, -1L)
})

# The class entries cover every link the package ships. They are compared on
# each component's own scale (expect_agrees_on_scale(), helper-links.R): the
# C side and the R methods run the same expressions, but where R evaluates a
# product and a sum separately a compiler may contract them into one
# multiply-add, so the identity of the two is a property of one compiler.
test_that("the class entries mirror the R methods of every link", {
  real_line <- c(-30, -5, -1, -0.2, 0.4, 2, 6, 30)
  positive <- c(0.01, 0.3, 1, 2.5, 40)
  for (nm in names(all_links())) {
    l <- all_links()[[nm]]
    rt <- link_scalar_route(l)
    expect_false(is.null(rt), label = nm)
    eta <- if (nm %in% c("inverse", "inverse_sq", "power_half", "power_2")) {
      positive
    } else real_line
    pr <- lf7_class_probe(rt$name, rt$par, eta)
    expect_gte(pr$id, 0L)
    cls <- S7::S7_class(l)
    ref <- list(h = S7::method(linkinv, cls)(l, eta),
                h1 = S7::method(dlinkinv, cls)(l, eta),
                h2 = S7::method(d2linkinv, cls)(l, eta))
    for (part in names(ref)) {
      if (all(ref[[part]] == 0)) {
        expect_identical(pr[[part]], as.numeric(ref[[part]]))
      } else {
        expect_agrees_on_scale(pr[[part]], as.numeric(ref[[part]]),
                               label = paste(nm, part))
      }
    }
  }
})

# The entry to the fourth order, lf7_inv1234p(), against linkinv() and the
# four derivative methods, on the grids above and at the extremes where
# rescue_pow() and the far branch of the logistic kernels take over.
test_that("the class entries to the fourth order mirror the R methods of every link", {
  real_line <- c(-800, -650, -30, -5, -1, -0.2, 0.4, 2, 6, 30, 650, 800)
  positive <- c(1e-200, 0.01, 0.3, 1, 2.5, 40, 1e200)
  gens <- list(linkinv, dlinkinv, d2linkinv, d3linkinv, d4linkinv)
  for (nm in names(all_links())) {
    l <- all_links()[[nm]]
    rt <- link_scalar_route(l)
    eta <- if (nm %in% c("inverse", "inverse_sq", "power_half", "power_2")) {
      positive
    } else real_line
    pr <- lf7_class_probe4(rt$name, rt$par, eta)
    expect_gte(pr$id, 0L)
    cls <- S7::S7_class(l)
    for (k in seq_along(gens)) {
      ref <- as.numeric(suppressWarnings(S7::method(gens[[k]], cls)(l, eta)))
      got <- pr$h[, k]
      lab <- paste(nm, "order", k - 1)
      # the extremes reach infinities: those must be the same values, and
      # the finite entries are compared on their own scale
      fin <- is.finite(ref)
      expect_identical(is.finite(got), fin, label = paste(lab, "finite"))
      expect_identical(got[!fin], ref[!fin], label = paste(lab, "non-finite"))
      if (all(ref[fin] == 0)) {
        expect_identical(got[fin], ref[fin], label = lab)
      } else {
        expect_agrees_on_scale(got[fin], ref[fin], label = lab)
      }
    }
    # a missing eta is missing at every order
    expect_true(all(is.na(lf7_class_probe4(rt$name, rt$par, NA_real_)$h)),
                label = paste(nm, "NA"))
  }
})

# The link itself, lf7_fun1p(), against linkfun() of every class, at the
# values that linkinv() returns on the grids above: the compiled inner fit of
# statmodels7 composes the two as inverse_link_derivs() does.
test_that("the class entry of the link itself mirrors linkfun() of every link", {
  real_line <- c(-800, -650, -30, -5, -1, -0.2, 0.4, 2, 6, 30, 650, 800)
  positive <- c(1e-200, 0.01, 0.3, 1, 2.5, 40, 1e200)
  for (nm in names(all_links())) {
    l <- all_links()[[nm]]
    rt <- link_scalar_route(l)
    eta <- if (nm %in% c("inverse", "inverse_sq", "power_half", "power_2")) {
      positive
    } else real_line
    th <- linkinv(l, eta)
    ref <- suppressWarnings(linkfun(l, th))
    got <- lf7_fun_probe(rt$name, rt$par, th)
    expect_identical(got, as.numeric(ref), label = nm)
  }
})

test_that("a link outside the package has no scalar route", {
  Mine <- S7::new_class("MyLogitLink", parent = LogitLink)
  expect_null(link_scalar_route(Mine(link_name = "logit",
                                     link_bounds = c(0, 1))))
  expect_identical(lf7_class_probe("no-such-link", numeric(0), 0)$id, -1L)
})

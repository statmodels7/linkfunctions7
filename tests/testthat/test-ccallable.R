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

test_that("a link outside the package has no scalar route", {
  Mine <- S7::new_class("MyLogitLink", parent = LogitLink)
  expect_null(link_scalar_route(Mine(link_name = "logit",
                                     link_bounds = c(0, 1))))
  expect_identical(lf7_class_probe("no-such-link", numeric(0), 0)$id, -1L)
})

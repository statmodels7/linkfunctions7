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

test_that("linkinv_derivs() returns the order-specific generics bit for bit", {
  gens <- list(dlinkinv, d2linkinv, d3linkinv, d4linkinv, d5linkinv, d6linkinv)
  for (nm in names(all_links())) {
    lk <- all_links()[[nm]]
    b <- eta_bounds(lk)
    lo <- if (is.finite(b[1])) b[1] + 0.1 else -3
    hi <- if (is.finite(b[2])) b[2] - 0.1 else 3
    eta <- c(seq(lo, hi, length.out = 23), NA)
    for (k in 1:6) {
      got <- linkinv_derivs(lk, eta, order = k)
      expect_length(got, k)
      for (j in seq_len(k)) {
        expect_identical(got[[j]], gens[[j]](lk, eta), info = paste(nm, k, j))
      }
    }
  }
})

test_that("the log link's method covers the floor and the extremes", {
  lk <- log_link()
  eta <- c(-1000, -800, -1e-300, 0, 1, 700, 710, Inf, -Inf, NaN, NA)
  got <- linkinv_derivs(lk, eta, order = 5)
  for (j in 1:5) {
    expect_identical(got[[j]], list(dlinkinv, d2linkinv, d3linkinv, d4linkinv,
                                    d5linkinv)[[j]](lk, eta))
  }
})

test_that("a subclass of the log link is read through its own methods", {
  Shifted <- S7::new_class("ShiftedLog", parent = linkfunctions7:::LogLink)
  S7::method(d2linkinv, Shifted) <- function(x, eta) 2 * exp(eta)
  lk <- Shifted(link_name = "shifted", link_bounds = c(0, Inf),
                link_params = list())
  got <- linkinv_derivs(lk, c(0, 1), order = 3)
  expect_identical(got[[2]], 2 * exp(c(0, 1)))
  expect_identical(got[[1]], dlinkinv(lk, c(0, 1)))
})

test_that("linkinv_derivs() rejects an order outside 1 to 6", {
  for (bad in list(0, 7, 2.5, NA, c(1, 2), "2")) {
    expect_error(linkinv_derivs(log_link(), 1, order = bad),
                 "must be one of 1 to 6")
  }
})

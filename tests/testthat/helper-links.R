# Every link the package ships, in one place. Lives in a helper rather than in
# a test file because more than one file needs it, and testthat gives each test
# file its own environment: a function defined at the top of one is not visible
# from another.

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

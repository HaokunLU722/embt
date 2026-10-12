test_that("covariate adjustment returns both versions and diagnostics", {
  dat <- make_embt_data()
  set.seed(88)
  fit <- embt(
    outcome ~ treatment, data = dat,
    adjust = ~ x1 + x2, B = 29
  )

  expect_named(fit$tests, c("unadjusted", "adjusted"))
  expect_equal(nrow(fit$adjustment$theta), 2L)
  expect_equal(ncol(fit$adjustment$theta), fit$grid_size)
  expect_lt(fit$adjustment$variance_identity_error, 1e-10)
  expect_true(all(is.finite(fit$process$adjusted_centered_ete)))
})

test_that("rescaling a fixed basis does not change adjusted results", {
  dat <- make_embt_data()
  W <- as.matrix(dat[c("x1", "x2")])
  set.seed(31)
  fit1 <- embt(dat$outcome, dat$treatment, adjust = W, B = 19)
  set.seed(31)
  fit2 <- embt(
    dat$outcome, dat$treatment,
    adjust = sweep(W, 2, c(100, 0.01), "*"), B = 19
  )

  expect_equal(
    fit1$process$adjusted_centered_ete,
    fit2$process$adjusted_centered_ete,
    tolerance = 1e-10
  )
  expect_equal(fit1$tests$adjusted$p_value,
               fit2$tests$adjusted$p_value)
})


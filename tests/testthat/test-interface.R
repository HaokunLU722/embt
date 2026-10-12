test_that("default and formula interfaces agree", {
  dat <- make_embt_data()

  set.seed(19)
  vector_fit <- embt(
    dat$outcome, dat$treatment,
    adjust = dat[c("x1", "x2")], B = 19
  )
  set.seed(19)
  formula_fit <- embt(
    outcome ~ treatment, data = dat,
    adjust = ~ x1 + x2, B = 19
  )

  expect_equal(vector_fit$process, formula_fit$process, tolerance = 1e-12)
  expect_equal(
    vector_fit$tests$adjusted$p_value,
    formula_fit$tests$adjusted$p_value
  )
  expect_equal(formula_fit$omitted, 0L)
})

test_that("factor treatment levels determine the comparison direction", {
  dat <- make_embt_data()
  dat$treatment <- factor(dat$treatment, levels = c(0, 1),
                          labels = c("Control", "Treatment"))
  set.seed(5)
  fit <- embt(outcome ~ treatment, data = dat, B = 19)

  expect_equal(fit$treatment_levels, c("Control", "Treatment"))
  expect_equal(fit$n0, fit$n1)
  expect_gt(fit$mean_shift, 0)
})

test_that("formula method handles complete cases consistently", {
  dat <- make_embt_data()
  dat$x1[3] <- NA_real_
  dat$outcome[7] <- NA_real_
  set.seed(8)
  fit <- embt(
    outcome ~ treatment, data = dat,
    adjust = ~ x1 + x2, B = 19
  )

  expect_equal(fit$omitted, 2L)
  expect_equal(fit$n, nrow(dat) - 2L)
})

test_that("invalid inputs fail clearly", {
  dat <- make_embt_data()
  expect_error(embt(dat$outcome, rep(1, nrow(dat)), B = 19), "coded as 0")
  expect_error(embt(dat$outcome, dat$treatment, tau = c(0.6, 0.9), B = 19),
               "contain 0.5")
  expect_error(embt(dat$outcome, dat$treatment, B = 5), "at least 19")
  expect_error(
    embt(dat$outcome, dat$treatment, adjust = rep(1, nrow(dat)), B = 19),
    "constant or numerically degenerate"
  )
})


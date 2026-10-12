test_that("the test is translation invariant", {
  dat <- make_embt_data()
  set.seed(101)
  original <- embt(dat$outcome, dat$treatment, B = 29)
  set.seed(101)
  shifted <- embt(dat$outcome + 1e4, dat$treatment, B = 29)

  expect_equal(
    original$process$centered_ete,
    shifted$process$centered_ete,
    tolerance = 1e-9
  )
  expect_equal(original$tests$unadjusted$p_value,
               shifted$tests$unadjusted$p_value)
})

test_that("an exact sample location shift gives a zero centered curve", {
  set.seed(4)
  y0 <- rnorm(40)
  y1 <- y0 + 2
  y <- c(y0, y1)
  d <- rep(0:1, each = 40)
  set.seed(11)
  fit <- embt(y, d, B = 19)

  expect_equal(fit$process$centered_ete, rep(0, fit$grid_size),
               tolerance = 1e-12)
  expect_equal(unname(fit$tests$unadjusted$statistic), rep(0, 4),
               tolerance = 1e-12)
  expect_equal(unname(fit$tests$unadjusted$p_value), rep(1, 4))
})

test_that("the exact regular-process factorization is numerically accurate", {
  dat <- make_embt_data()
  set.seed(3)
  fit <- embt(dat$outcome, dat$treatment, B = 19)
  expect_lt(fit$diagnostics$factorization_error, 1e-10)
  expect_equal(fit$process$centered_ete[fit$tau_grid == 0.5], 0)
})

test_that("Mammen multipliers have the required analytical moments", {
  parameters <- embt:::.embt_mammen_parameters()
  moments <- vapply(1:3, function(power) {
    sum(parameters$probabilities * parameters$values^power)
  }, numeric(1))
  expect_equal(moments, c(0, 1, 1), tolerance = 1e-14)
})


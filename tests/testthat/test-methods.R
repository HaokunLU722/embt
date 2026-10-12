test_that("summary, confidence-band, and plotting methods work", {
  dat <- make_embt_data()
  set.seed(22)
  fit <- embt(
    outcome ~ treatment, data = dat,
    adjust = ~ x1, B = 19
  )

  summary_fit <- summary(fit)
  expect_s3_class(summary_fit, "summary.embt")
  expect_equal(nrow(summary_fit$coefficients), 8L)

  band <- confint(fit)
  adjusted_band <- confint(fit, parm = "adjusted")
  expect_equal(nrow(band), fit$grid_size)
  expect_true(all(band$lower <= band$estimate))
  expect_true(all(adjusted_band$upper >= adjusted_band$estimate))

  filename <- tempfile(fileext = ".pdf")
  grDevices::pdf(filename)
  expect_invisible(plot(fit, version = "both"))
  grDevices::dev.off()
  expect_true(file.exists(filename))
})


.embt_bootstrap_process <- function(fit, V1, V0, chi1, chi0) {
  sqrt(1 - fit$lambda) * (V1 %*% chi1) / sqrt(fit$n1) -
    sqrt(fit$lambda) * (V0 %*% chi0) / sqrt(fit$n0)
}

.embt_observed_statistics <- function(S, U, tau_grid) {
  raw <- .embt_process_statistics(S, tau_grid)
  studentized <- .embt_process_statistics(U, tau_grid)
  c(
    KS = unname(raw["KS"]),
    CvM = unname(raw["CvM"]),
    S_KS = unname(studentized["KS"]),
    S_CvM = unname(studentized["CvM"])
  )
}

.embt_bootstrap_statistics <- function(A_star, c_tau, variance, tau_grid) {
  raw <- .embt_process_statistics(sweep(A_star, 2L, c_tau, "*"), tau_grid)
  studentized <- .embt_process_statistics(
    sweep(A_star, 2L, sqrt(variance), "/"),
    tau_grid
  )
  cbind(
    KS = raw[, "KS"],
    CvM = raw[, "CvM"],
    S_KS = studentized[, "KS"],
    S_CvM = studentized[, "CvM"]
  )
}

.embt_p_values <- function(observed, bootstrap) {
  vapply(names(observed), function(name) {
    (1 + sum(bootstrap[, name] >= observed[name])) / (nrow(bootstrap) + 1)
  }, numeric(1))
}

.embt_run_bootstrap <- function(fit, adjustment, B, control) {
  unadjusted <- matrix(
    NA_real_, B, 4L,
    dimnames = list(NULL, c("KS", "CvM", "S_KS", "S_CvM"))
  )
  adjusted <- if (is.null(adjustment)) NULL else unadjusted
  starts <- seq.int(1L, B, by = control$batch_size)

  for (start in starts) {
    end <- min(B, start + control$batch_size - 1L)
    current_B <- end - start + 1L
    V1 <- .embt_draw_mammen_matrix(current_B, fit$n1)
    V0 <- .embt_draw_mammen_matrix(current_B, fit$n0)
    A_star <- .embt_bootstrap_process(
      fit, V1, V0, fit$group1$chi, fit$group0$chi
    )
    unadjusted[start:end, ] <- .embt_bootstrap_statistics(
      A_star, fit$c_tau, fit$v2, fit$tau_grid
    )

    if (!is.null(adjustment)) {
      A_star_adjusted <- .embt_bootstrap_process(
        fit, V1, V0, adjustment$chi1, adjustment$chi0
      )
      adjusted[start:end, ] <- .embt_bootstrap_statistics(
        A_star_adjusted, fit$c_tau, adjustment$v2, fit$tau_grid
      )
    }

    if (control$verbose) {
      message(sprintf("Completed %d of %d multiplier draws.", end, B))
    }
  }

  list(unadjusted = unadjusted, adjusted = adjusted)
}

.embt_test_block <- function(S, U, bootstrap, tau_grid) {
  observed <- .embt_observed_statistics(S, U, tau_grid)
  list(
    statistic = observed,
    p_value = .embt_p_values(observed, bootstrap),
    bootstrap = bootstrap
  )
}

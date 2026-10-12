.embt_expectile_sorted <- function(y_sorted, cumulative_sum, tau) {
  if (tau == 0.5) {
    return(mean(y_sorted))
  }

  n <- length(y_sorted)
  k <- 0:n
  total <- cumulative_sum[n + 1L]
  numerator <- (1 - tau) * cumulative_sum +
    tau * (total - cumulative_sum)
  denominator <- (1 - tau) * k + tau * (n - k)
  candidates <- numerator / denominator
  lower <- c(-Inf, y_sorted)
  upper <- c(y_sorted, Inf)
  tolerance <- 64 * .Machine$double.eps * max(1, max(abs(y_sorted)))
  admissible <- which(
    candidates >= lower - tolerance & candidates <= upper + tolerance
  )

  snap_to_observation <- function(value) {
    nearest <- y_sorted[which.min(abs(y_sorted - value))]
    if (abs(nearest - value) <= tolerance) nearest else value
  }

  if (length(admissible) == 1L) {
    return(snap_to_observation(candidates[admissible]))
  }
  if (length(admissible) > 1L) {
    score <- vapply(candidates[admissible], function(value) {
      abs(sum(2 * abs(tau - (y_sorted <= value)) * (y_sorted - value)))
    }, numeric(1))
    return(snap_to_observation(candidates[admissible[which.min(score)]]))
  }

  score <- function(value) {
    sum(2 * abs(tau - (y_sorted <= value)) * (y_sorted - value))
  }
  stats::uniroot(
    score,
    interval = range(y_sorted),
    tol = .Machine$double.eps^0.75
  )$root
}

.embt_sample_expectiles <- function(y, tau_grid) {
  y_sorted <- sort(y)
  cumulative_sum <- c(0, cumsum(y_sorted))
  values <- vapply(tau_grid, function(tau) {
    .embt_expectile_sorted(y_sorted, cumulative_sum, tau)
  }, numeric(1))
  values[tau_grid == 0.5] <- mean(y)
  values
}

.embt_regular_component <- function(y, tau_grid) {
  n <- length(y)
  y_mean <- mean(y)
  y_sorted <- sort(y)
  expectiles <- .embt_sample_expectiles(y, tau_grid)
  distribution <- findInterval(expectiles, y_sorted) / n
  sensitivity <- 2 * (tau_grid + (1 - 2 * tau_grid) * distribution)
  m <- numeric(length(tau_grid))
  chi <- matrix(NA_real_, nrow = n, ncol = length(tau_grid))

  for (j in seq_along(tau_grid)) {
    absolute_residual <- abs(y - expectiles[j])
    m[j] <- mean(absolute_residual)
    a <- 1 - 2 * distribution[j]
    chi[, j] <- (
      absolute_residual - m[j] - a * (y - y_mean)
    ) / sensitivity[j]
  }

  colnames(chi) <- format(tau_grid, trim = TRUE)
  list(
    n = n,
    mean = y_mean,
    expectiles = expectiles,
    distribution = distribution,
    sensitivity = sensitivity,
    m = m,
    chi = chi
  )
}

.embt_regular_fit <- function(y1, y0, tau_grid) {
  group1 <- .embt_regular_component(y1, tau_grid)
  group0 <- .embt_regular_component(y0, tau_grid)
  n1 <- length(y1)
  n0 <- length(y0)
  n <- n1 + n0
  lambda <- n1 / n
  effective_n <- n1 * n0 / n
  c_tau <- 2 * tau_grid - 1
  delta_hat <- mean(y1) - mean(y0)
  centered_ete <- group1$expectiles - group0$expectiles - delta_hat
  centered_ete[tau_grid == 0.5] <- 0

  A <- sqrt(effective_n) * (group1$m - group0$m)
  S <- c_tau * A
  factorization_error <- max(abs(S / sqrt(effective_n) - centered_ete))

  v2 <- (1 - lambda) * colMeans(group1$chi^2) +
    lambda * colMeans(group0$chi^2)
  variance_floor <- 64 * .Machine$double.eps * max(1, max(v2))
  .embt_assert(
    all(is.finite(v2) & v2 > variance_floor),
    "The regular-process variance is not positive over the expectile range."
  )

  list(
    tau_grid = tau_grid,
    c_tau = c_tau,
    n1 = n1,
    n0 = n0,
    n = n,
    lambda = lambda,
    effective_n = effective_n,
    delta_hat = delta_hat,
    group1 = group1,
    group0 = group0,
    centered_ete = centered_ete,
    A = A,
    S = S,
    U = A / sqrt(v2),
    v2 = v2,
    factorization_error = factorization_error
  )
}


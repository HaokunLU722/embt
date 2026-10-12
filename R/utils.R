.embt_assert <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
  invisible(TRUE)
}

.embt_validate_outcome <- function(y) {
  .embt_assert(is.numeric(y), "`outcome` must be numeric.")
  .embt_assert(length(y) >= 4L, "`outcome` must contain at least four observations.")
  .embt_assert(all(is.finite(y)), "`outcome` must contain only finite values.")
  as.numeric(y)
}

.embt_binary_treatment <- function(treatment) {
  .embt_assert(length(treatment) >= 4L, "`treatment` is too short.")
  .embt_assert(!anyNA(treatment), "`treatment` must not contain missing values.")

  if (is.logical(treatment)) {
    z <- as.integer(treatment)
    labels <- c("FALSE", "TRUE")
  } else if (is.factor(treatment) || is.character(treatment)) {
    treatment <- droplevels(factor(treatment))
    .embt_assert(
      nlevels(treatment) == 2L,
      "A factor or character `treatment` must have exactly two levels."
    )
    z <- as.integer(treatment) - 1L
    labels <- levels(treatment)
  } else if (is.numeric(treatment) || is.integer(treatment)) {
    values <- sort(unique(treatment))
    .embt_assert(
      identical(as.numeric(values), c(0, 1)),
      "A numeric `treatment` must be coded as 0 for control and 1 for treatment."
    )
    z <- as.integer(treatment)
    labels <- c("0", "1")
  } else {
    stop("`treatment` must be logical, binary numeric, character, or a factor.",
         call. = FALSE)
  }

  .embt_assert(all(tabulate(z + 1L, nbins = 2L) >= 2L),
               "Each treatment arm must contain at least two observations.")
  list(z = z, labels = labels)
}

.embt_tau_grid <- function(tau, grid_size) {
  .embt_assert(is.numeric(tau) && length(tau) == 2L && all(is.finite(tau)),
               "`tau` must be a finite numeric vector of length two.")
  tau <- sort(as.numeric(tau))
  .embt_assert(tau[1L] > 0 && tau[2L] < 1,
               "The endpoints of `tau` must lie strictly between zero and one.")
  .embt_assert(tau[1L] < 0.5 && tau[2L] > 0.5,
               "The expectile range must contain 0.5 in its interior.")

  grid <- seq(tau[1L], tau[2L], length.out = grid_size)
  nearest <- which.min(abs(grid - 0.5))
  grid[nearest] <- 0.5
  grid <- sort(grid)
  .embt_assert(!is.unsorted(grid, strictly = TRUE),
               "The numerical expectile grid is not strictly increasing.")
  grid
}

.embt_cvm_weights <- function(tau_grid) {
  delta <- diff(tau_grid)
  weights <- numeric(length(tau_grid))
  weights[1L] <- delta[1L] / 2
  weights[length(weights)] <- delta[length(delta)] / 2
  if (length(weights) > 2L) {
    weights[2L:(length(weights) - 1L)] <-
      (delta[-length(delta)] + delta[-1L]) / 2
  }
  weights / diff(range(tau_grid))
}

.embt_process_statistics <- function(process, tau_grid) {
  weights <- .embt_cvm_weights(tau_grid)
  if (is.null(dim(process))) {
    return(c(KS = max(abs(process)), CvM = sum(weights * process^2)))
  }
  cbind(
    KS = apply(abs(process), 1L, max),
    CvM = drop((process^2) %*% weights)
  )
}

.embt_quantile <- function(x, probability) {
  unname(stats::quantile(x, probs = probability, type = 1, names = FALSE))
}

.embt_mammen_parameters <- function() {
  kappa <- (sqrt(5) + 1) / 2
  list(
    values = c(1 - kappa, kappa),
    probabilities = c(kappa / sqrt(5), 1 - kappa / sqrt(5))
  )
}

.embt_draw_mammen_matrix <- function(B, n) {
  parameters <- .embt_mammen_parameters()
  uniforms <- stats::runif(B * n)
  values <- ifelse(
    uniforms < parameters$probabilities[1L],
    parameters$values[1L],
    parameters$values[2L]
  )
  matrix(values, nrow = B, ncol = n)
}

.embt_match_version <- function(version, has_adjustment) {
  if (length(version) > 1L) {
    version <- version[1L]
  }
  choices <- if (has_adjustment) c("adjusted", "unadjusted") else "unadjusted"
  match.arg(version, choices)
}

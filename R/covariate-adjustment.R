.embt_prepare_basis <- function(W, treatment) {
  if (is.data.frame(W)) {
    W <- stats::model.matrix(~ . - 1, data = W)
  }
  if (is.null(dim(W))) {
    W <- matrix(W, ncol = 1L)
  }
  W <- as.matrix(W)
  storage.mode(W) <- "double"
  .embt_assert(nrow(W) == length(treatment),
               "The adjustment basis and outcome must have the same number of rows.")
  .embt_assert(ncol(W) >= 1L,
               "The adjustment basis must contain at least one non-intercept column.")
  .embt_assert(all(is.finite(W)),
               "The adjustment basis must contain only finite values.")

  if (is.null(colnames(W))) {
    colnames(W) <- paste0("basis", seq_len(ncol(W)))
  }

  scales <- apply(W, 2L, stats::sd)
  tolerance <- sqrt(.Machine$double.eps) * max(1, max(abs(W)))
  keep_scale <- is.finite(scales) & scales > tolerance
  dropped <- colnames(W)[!keep_scale]
  W <- W[, keep_scale, drop = FALSE]
  scales <- scales[keep_scale]
  .embt_assert(ncol(W) >= 1L,
               "All adjustment-basis columns are constant or numerically degenerate.")
  W <- sweep(W, 2L, scales, "/")

  W1 <- W[treatment == 1L, , drop = FALSE]
  W0 <- W[treatment == 0L, , drop = FALSE]
  W1_centered <- sweep(W1, 2L, colMeans(W1), "-")
  W0_centered <- sweep(W0, 2L, colMeans(W0), "-")
  centered <- rbind(W1_centered, W0_centered)
  decomposition <- qr(centered, tol = 1e-10)
  if (decomposition$rank < ncol(W)) {
    .embt_assert(decomposition$rank >= 1L,
                 "The adjustment basis has no within-arm variation.")
    keep_rank <- sort(decomposition$pivot[seq_len(decomposition$rank)])
    dropped <- c(dropped, colnames(W)[-keep_rank])
    W <- W[, keep_rank, drop = FALSE]
    warning(
      "Redundant adjustment-basis columns were removed: ",
      paste(unique(dropped), collapse = ", "),
      call. = FALSE
    )
  } else if (length(dropped)) {
    warning(
      "Constant adjustment-basis columns were removed: ",
      paste(unique(dropped), collapse = ", "),
      call. = FALSE
    )
  }

  list(W = W, names = colnames(W), dropped = unique(dropped))
}

.embt_adjustment_fit <- function(fit, W, treatment) {
  basis <- .embt_prepare_basis(W, treatment)
  W1 <- basis$W[treatment == 1L, , drop = FALSE]
  W0 <- basis$W[treatment == 0L, , drop = FALSE]
  W1_centered <- sweep(W1, 2L, colMeans(W1), "-")
  W0_centered <- sweep(W0, 2L, colMeans(W0), "-")

  sigma <- (1 - fit$lambda) * crossprod(W1_centered) / fit$n1 +
    fit$lambda * crossprod(W0_centered) / fit$n0
  condition_number <- kappa(sigma, exact = TRUE)
  .embt_assert(is.finite(condition_number),
               "The adjustment covariance matrix is singular.")

  C <- (1 - fit$lambda) *
    crossprod(W1_centered, fit$group1$chi) / fit$n1 +
    fit$lambda * crossprod(W0_centered, fit$group0$chi) / fit$n0
  theta <- solve(sigma, C)
  imbalance <- sqrt(fit$effective_n) * (colMeans(W1) - colMeans(W0))
  A <- fit$A - drop(imbalance %*% theta)
  S <- fit$c_tau * A

  chi1 <- fit$group1$chi - W1_centered %*% theta
  chi0 <- fit$group0$chi - W0_centered %*% theta
  v2 <- (1 - fit$lambda) * colMeans(chi1^2) +
    fit$lambda * colMeans(chi0^2)
  variance_floor <- 64 * .Machine$double.eps * max(1, max(v2))
  .embt_assert(
    all(is.finite(v2) & v2 > variance_floor),
    "The adjusted regular-process variance is not positive over the expectile range."
  )

  explained <- colSums(C * theta)
  list(
    basis_names = basis$names,
    dropped_basis = basis$dropped,
    W1_centered = W1_centered,
    W0_centered = W0_centered,
    sigma = sigma,
    C = C,
    theta = theta,
    imbalance = imbalance,
    condition_number = condition_number,
    A = A,
    S = S,
    U = A / sqrt(v2),
    v2 = v2,
    chi1 = chi1,
    chi0 = chi0,
    centered_ete = S / sqrt(fit$effective_n),
    variance_reduction = 1 - v2 / fit$v2,
    projection_R2 = explained / fit$v2,
    variance_identity_error = max(abs(v2 - (fit$v2 - explained)))
  )
}

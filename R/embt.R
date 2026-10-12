#' Expectile multiplier-bootstrap tests of treatment-effect heterogeneity
#'
#' Tests whether treatment and control outcome distributions differ only by a
#' location shift. The function jointly computes raw and pointwise-studentized
#' Kolmogorov--Smirnov and Cramer--von Mises tests. If a fixed covariate basis is
#' supplied, covariate-adjusted versions are computed from the same multiplier
#' draws.
#'
#' @param formula A two-sided formula of the form `outcome ~ treatment`, or a
#'   numeric outcome vector when using the default method.
#' @param data A data frame used with the formula method.
#' @param adjust An optional one-sided formula specifying a fixed,
#'   finite-dimensional covariate basis, such as `~ x1 + x2`. With the default
#'   method, supply a numeric vector, matrix, or data frame.
#' @param tau A length-two vector giving the lower and upper expectile indices.
#'   The interval must contain 0.5 in its interior.
#' @param B Number of Mammen multiplier-bootstrap draws.
#' @param control Advanced numerical controls created by [embt_control()].
#' @param ... Arguments passed to methods.
#'
#' @details A numeric treatment must be coded zero for control and one for
#'   treatment. For a two-level factor or character treatment, the first level
#'   is control and the second is treatment. The default implementation uses 81
#'   numerical expectile coordinates and does not estimate densities, sparsity
#'   functions, or score functions.
#'
#' @return An object of class `embt` containing process estimates, four tests,
#'   multiplier statistics, and, when requested, covariate-adjustment
#'   diagnostics.
#'
#' @examples
#' set.seed(42)
#' n <- 80
#' d <- rep(0:1, each = n / 2)
#' x <- rnorm(n)
#' y <- 1 + d + (1 + 0.3 * d) * rnorm(n) + 0.4 * x
#'
#' fit <- embt(y, d, B = 49)
#' summary(fit)
#'
#' fit_ca <- embt(y, d, adjust = x, B = 49)
#' summary(fit_ca)
#'
#' @export
embt <- function(formula, ...) {
  UseMethod("embt")
}

#' @rdname embt
#' @export
embt.formula <- function(formula, data, adjust = NULL,
                         tau = c(0.10, 0.90), B = 999L,
                         control = embt_control(), ...) {
  .embt_assert(inherits(formula, "formula"), "`formula` must be a formula.")
  .embt_assert(length(formula) == 3L,
               "`formula` must be two-sided, for example `outcome ~ treatment`.")
  .embt_assert(!missing(data) && is.data.frame(data),
               "`data` must be supplied as a data frame.")

  outcome_frame <- stats::model.frame(
    formula, data = data, na.action = stats::na.pass
  )
  .embt_assert(ncol(outcome_frame) == 2L,
               "The right-hand side of `formula` must specify one treatment variable.")

  adjustment_frame <- NULL
  if (!is.null(adjust)) {
    .embt_assert(inherits(adjust, "formula") && length(adjust) == 2L,
                 "`adjust` must be a one-sided formula, for example `~ x1 + x2`.")
    adjustment_frame <- stats::model.frame(
      adjust, data = data, na.action = stats::na.pass
    )
  }

  complete <- stats::complete.cases(outcome_frame)
  if (!is.null(adjustment_frame)) {
    complete <- complete & stats::complete.cases(adjustment_frame)
  }
  .embt_assert(any(complete), "No complete observations remain.")
  omitted <- sum(!complete)

  y <- stats::model.response(outcome_frame)[complete]
  treatment <- outcome_frame[[2L]][complete]
  W <- NULL
  if (!is.null(adjust)) {
    W <- stats::model.matrix(adjust, data = data[complete, , drop = FALSE])
    intercept <- which(colnames(W) == "(Intercept)")
    if (length(intercept)) {
      W <- W[, -intercept, drop = FALSE]
    }
    .embt_assert(ncol(W) >= 1L,
                 "`adjust` must contain at least one non-intercept basis term.")
  }

  object <- embt.default(
    y,
    treatment = treatment,
    adjust = W,
    tau = tau,
    B = B,
    control = control
  )
  object$call <- match.call()
  object$formula <- formula
  object$adjust_formula <- adjust
  object$omitted <- omitted
  object
}

#' @rdname embt
#' @param treatment A binary treatment vector for the default method.
#' @export
embt.default <- function(formula, treatment, adjust = NULL,
                         tau = c(0.10, 0.90), B = 999L,
                         control = embt_control(), ...) {
  y <- .embt_validate_outcome(formula)
  .embt_assert(length(treatment) == length(y),
               "`outcome` and `treatment` must have the same length.")
  treatment_info <- .embt_binary_treatment(treatment)
  z <- treatment_info$z
  control <- .embt_validate_control(control)

  B <- as.integer(B)
  .embt_assert(length(B) == 1L && is.finite(B) && B >= 19L,
               "`B` must be an integer of at least 19.")
  tau_grid <- .embt_tau_grid(tau, control$grid_size)
  y1 <- y[z == 1L]
  y0 <- y[z == 0L]
  fit <- .embt_regular_fit(y1, y0, tau_grid)

  adjustment <- if (is.null(adjust)) {
    NULL
  } else {
    .embt_adjustment_fit(fit, adjust, z)
  }
  bootstrap <- .embt_run_bootstrap(fit, adjustment, B, control)
  tests <- list(
    unadjusted = .embt_test_block(
      fit$S, fit$U, bootstrap$unadjusted, tau_grid
    )
  )
  if (!is.null(adjustment)) {
    tests$adjusted <- .embt_test_block(
      adjustment$S, adjustment$U, bootstrap$adjusted, tau_grid
    )
  }

  process <- data.frame(
    tau = tau_grid,
    centered_ete = fit$centered_ete,
    regular_process = fit$A,
    regular_variance = fit$v2,
    studentized_process = fit$U
  )
  if (!is.null(adjustment)) {
    process$adjusted_centered_ete <- adjustment$centered_ete
    process$adjusted_regular_process <- adjustment$A
    process$adjusted_regular_variance <- adjustment$v2
    process$adjusted_studentized_process <- adjustment$U
    process$variance_reduction <- adjustment$variance_reduction
  }

  adjustment_output <- NULL
  if (!is.null(adjustment)) {
    adjustment_output <- adjustment[c(
      "basis_names", "dropped_basis", "theta", "sigma", "C",
      "condition_number", "variance_reduction", "projection_R2",
      "variance_identity_error"
    )]
  }

  structure(
    list(
      call = match.call(),
      tau = range(tau_grid),
      tau_grid = tau_grid,
      grid_size = length(tau_grid),
      B = B,
      n = fit$n,
      n1 = fit$n1,
      n0 = fit$n0,
      effective_n = fit$effective_n,
      treatment_levels = treatment_info$labels,
      mean_shift = fit$delta_hat,
      process = process,
      tests = tests,
      adjustment = adjustment_output,
      diagnostics = list(
        factorization_error = fit$factorization_error,
        multiplier = "Mammen two-point"
      ),
      omitted = 0L
    ),
    class = "embt"
  )
}


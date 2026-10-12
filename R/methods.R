.embt_summary_table <- function(object, level) {
  rows <- lapply(names(object$tests), function(version) {
    test <- object$tests[[version]]
    data.frame(
      version = version,
      test = c("KS", "CvM", "S-KS", "S-CvM"),
      statistic = unname(test$statistic),
      critical_value = vapply(seq_along(test$statistic), function(j) {
        .embt_quantile(test$bootstrap[, j], level)
      }, numeric(1)),
      p_value = unname(test$p_value),
      row.names = NULL
    )
  })
  do.call(rbind, rows)
}

#' @export
print.embt <- function(x, ...) {
  cat("Expectile multiplier-bootstrap test\n")
  cat(sprintf(
    "Treatment: %s (n = %d); control: %s (n = %d)\n",
    x$treatment_levels[2L], x$n1, x$treatment_levels[1L], x$n0
  ))
  cat(sprintf(
    "Expectile range: [%.3f, %.3f]; grid: %d; Mammen draws: %d\n",
    x$tau[1L], x$tau[2L], x$grid_size, x$B
  ))
  cat(sprintf("Estimated mean shift: %.6g\n", x$mean_shift))
  print(.embt_summary_table(x, 0.95)[, c("version", "test", "p_value")],
        row.names = FALSE, digits = 4)
  invisible(x)
}

#' Summarize an eMBT fit
#'
#' @param object An object returned by [embt()].
#' @param level Confidence level used for multiplier critical values.
#' @param ... Unused.
#'
#' @return An object of class `summary.embt`.
#' @export
summary.embt <- function(object, level = 0.95, ...) {
  .embt_assert(is.numeric(level) && length(level) == 1L &&
                 is.finite(level) && level > 0 && level < 1,
               "`level` must lie strictly between zero and one.")
  structure(
    list(
      call = object$call,
      n = object$n,
      n1 = object$n1,
      n0 = object$n0,
      treatment_levels = object$treatment_levels,
      mean_shift = object$mean_shift,
      tau = object$tau,
      grid_size = object$grid_size,
      B = object$B,
      level = level,
      coefficients = .embt_summary_table(object, level),
      adjustment = object$adjustment,
      omitted = object$omitted
    ),
    class = "summary.embt"
  )
}

#' @export
print.summary.embt <- function(x, ...) {
  cat("Expectile multiplier-bootstrap tests\n\n")
  cat(sprintf(
    "Treatment: %s (n = %d); control: %s (n = %d)\n",
    x$treatment_levels[2L], x$n1, x$treatment_levels[1L], x$n0
  ))
  cat(sprintf(
    "Expectile range: [%.3f, %.3f]; numerical grid: %d points\n",
    x$tau[1L], x$tau[2L], x$grid_size
  ))
  cat(sprintf("Mammen multiplier draws: %d\n", x$B))
  cat(sprintf("Estimated mean shift: %.6g\n", x$mean_shift))
  if (x$omitted > 0L) {
    cat(sprintf("Observations omitted because of missing values: %d\n", x$omitted))
  }
  if (!is.null(x$adjustment)) {
    cat(sprintf(
      "Adjustment basis: %s\n",
      paste(x$adjustment$basis_names, collapse = ", ")
    ))
  }
  cat("\n")
  table <- x$coefficients
  table$statistic <- formatC(table$statistic, digits = 4, format = "fg")
  table$critical_value <- formatC(
    table$critical_value, digits = 4, format = "fg"
  )
  table$p_value <- format.pval(table$p_value, digits = 4, eps = 1 / (x$B + 1))
  names(table) <- c(
    "Version", "Test", "Statistic",
    sprintf("%.1f%% critical value", 100 * x$level), "p-value"
  )
  print(table, row.names = FALSE, right = FALSE)
  invisible(x)
}

#' Simultaneous confidence band for the centered ETE curve
#'
#' Constructs a raw KS-based simultaneous confidence band for the population
#' centered expectile treatment-effect curve.
#'
#' @param object An object returned by [embt()].
#' @param parm Either `"unadjusted"` or `"adjusted"`.
#' @param level Confidence level.
#' @param ... Unused.
#'
#' @return A data frame containing the expectile index, point estimate, and
#'   simultaneous lower and upper bounds.
#' @export
confint.embt <- function(object, parm = c("unadjusted", "adjusted"),
                         level = 0.95, ...) {
  parm <- .embt_match_version(parm, !is.null(object$adjustment))
  .embt_assert(is.numeric(level) && length(level) == 1L &&
                 is.finite(level) && level > 0 && level < 1,
               "`level` must lie strictly between zero and one.")
  test <- object$tests[[parm]]
  critical <- .embt_quantile(test$bootstrap[, "KS"], level)
  estimate <- if (parm == "adjusted") {
    object$process$adjusted_centered_ete
  } else {
    object$process$centered_ete
  }
  half_width <- critical / sqrt(object$effective_n)
  data.frame(
    tau = object$tau_grid,
    estimate = estimate,
    lower = estimate - half_width,
    upper = estimate + half_width,
    version = parm,
    level = level
  )
}

#' Plot an estimated centered ETE curve
#'
#' @param x An object returned by [embt()].
#' @param version Curve to display: `"unadjusted"`, `"adjusted"`, or `"both"`.
#' @param band Whether to add a raw KS simultaneous confidence band. Bands are
#'   available when a single version is plotted.
#' @param level Confidence level for the optional band.
#' @param ... Additional arguments passed to [graphics::plot()].
#'
#' @return The input object, invisibly.
#' @export
plot.embt <- function(x, version = c("unadjusted", "adjusted", "both"),
                      band = FALSE, level = 0.95, ...) {
  available_adjustment <- !is.null(x$adjustment)
  version <- match.arg(version)
  if (!available_adjustment && version != "unadjusted") {
    stop("Adjusted results are not available in this fit.", call. = FALSE)
  }
  if (isTRUE(band) && version == "both") {
    stop("Set `band = FALSE` when plotting both curves.", call. = FALSE)
  }

  tau <- x$tau_grid
  unadjusted <- x$process$centered_ete
  adjusted <- x$process$adjusted_centered_ete
  primary <- if (version == "adjusted") adjusted else unadjusted
  y_values <- primary
  confidence <- NULL
  if (isTRUE(band)) {
    confidence <- stats::confint(x, parm = version, level = level)
    y_values <- c(y_values, confidence$lower, confidence$upper)
  }
  if (version == "both") {
    y_values <- c(unadjusted, adjusted)
  }

  graphics::plot(
    tau, primary,
    type = "l", lwd = 2, col = "#B51F32",
    xlab = expression(tau),
    ylab = "Centered expectile treatment effect",
    ylim = range(y_values, finite = TRUE), ...
  )
  graphics::abline(h = 0, col = "grey65", lwd = 0.8)
  graphics::abline(v = 0.5, col = "grey75", lwd = 0.8)

  if (!is.null(confidence)) {
    graphics::lines(tau, confidence$lower, lty = 2, col = "grey35")
    graphics::lines(tau, confidence$upper, lty = 2, col = "grey35")
  }
  if (version == "both") {
    graphics::lines(tau, adjusted, lwd = 2, lty = 2, col = "#287DB5")
    graphics::legend(
      "topright",
      legend = c("Unadjusted", "Covariate adjusted"),
      col = c("#B51F32", "#287DB5"),
      lty = c(1, 2), lwd = 2, bty = "n"
    )
  }
  invisible(x)
}

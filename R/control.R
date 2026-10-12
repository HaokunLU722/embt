#' Advanced controls for eMBT
#'
#' These settings govern numerical implementation rather than statistical
#' smoothing or nuisance estimation. Most users should keep the defaults.
#'
#' @param grid_size Number of numerical expectile coordinates used over the
#'   requested range. The default is 81.
#' @param batch_size Number of multiplier draws processed in one batch.
#' @param verbose Whether to report multiplier-bootstrap progress.
#'
#' @return An object of class `embt_control`.
#' @export
embt_control <- function(grid_size = 81L, batch_size = 100L, verbose = FALSE) {
  grid_size <- as.integer(grid_size)
  batch_size <- as.integer(batch_size)
  .embt_assert(is.finite(grid_size) && grid_size >= 11L,
               "`grid_size` must be an integer of at least 11.")
  .embt_assert(is.finite(batch_size) && batch_size >= 1L,
               "`batch_size` must be a positive integer.")
  .embt_assert(is.logical(verbose) && length(verbose) == 1L && !is.na(verbose),
               "`verbose` must be TRUE or FALSE.")
  structure(
    list(grid_size = grid_size, batch_size = batch_size, verbose = verbose),
    class = "embt_control"
  )
}

.embt_validate_control <- function(control) {
  .embt_assert(inherits(control, "embt_control"),
               "`control` must be created by `embt_control()`.")
  control
}


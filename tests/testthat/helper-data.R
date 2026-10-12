make_embt_data <- function(n = 80L, heterogeneous = TRUE) {
  set.seed(7301)
  treatment <- rep(0:1, each = n / 2)
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.4)
  scale <- if (heterogeneous) 1 + 0.35 * treatment else 1
  outcome <- 0.5 + treatment + 0.4 * x1 + scale * rnorm(n)
  data.frame(outcome, treatment, x1, x2)
}


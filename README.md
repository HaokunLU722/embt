
# embt

<!-- badges: start -->
<!-- badges: end -->

`embt` implements expectile multiplier-bootstrap tests of whether
treatment and control outcome distributions differ only by a location
shift. The method does not estimate densities, sparsity functions, or
score functions.

## Installation

The development version can be installed directly from GitHub:

``` r
remotes::install_github("HaokunLU722/embt")
```

## Basic use

``` r
library(embt)

set.seed(123)
fit <- embt(outcome ~ treatment, data = experiment, B = 999)
summary(fit)
plot(fit)
```

One call returns the raw KS and CvM tests and their
pointwise-studentized counterparts. The default expectile range is
`[0.10, 0.90]`, evaluated on an internally fixed 81-point numerical
grid, and calibration uses Mammen two-point multipliers.

## Covariate adjustment

Provide a prespecified, finite-dimensional baseline basis with a
one-sided formula:

``` r
set.seed(123)
fit_ca <- embt(
  outcome ~ treatment,
  data = experiment,
  adjust = ~ baseline_outcome + age + factor(region),
  B = 999
)

summary(fit_ca)
plot(fit_ca, version = "both")
```

The adjusted and unadjusted tests use the same multiplier draws. No
smoothing or nonparametric nuisance estimation is introduced by the
adjustment.

## Simultaneous confidence band

``` r
band <- confint(fit, level = 0.95)
plot(fit, band = TRUE)
```

The confidence band is based on the raw KS multiplier statistic.

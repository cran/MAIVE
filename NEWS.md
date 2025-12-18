# MAIVE 0.1.11

*Released: 2025-12-18*

## Other Changes

* Disable instrumentation when Ns has no variation; avoid aliased-slope vcov indexing
* Guard first-stage F-test against rank-deficient vcovCR; add regression test

---


# MAIVE 0.1.10

*Released: 2025-12-02*

## Internal

* Update outdated cran files

---


# MAIVE 0.1.9

*Released: 2025-12-02*

## Bug Fixes

* Cran submission issues

---


# MAIVE 0.1.8

*Released: 2025-11-27*

## Bug Fixes

* Further issues in the ar calculation
* Ar SE usage


## Internal

* Fix outdated tests

---


# MAIVE 0.1.7

*Released: 2025-11-27*

## Internal

* Keep the ar tests more neutral
* Ar calculation - avoid banana projection

---


# MAIVE 0.1.6

*Released: 2025-11-26*

## New Features

* Add automatic news updates


## Bug Fixes

* Automatic news release


## Documentation

* Update release instructions docs

---


# MAIVE 0.0.4

## Initial CRAN submission

* Implemented MAIVE (Meta-Analysis Instrumental Variable Estimator) for addressing spurious precision in meta-analysis
* Core functions:
  * `maive()`: Main function implementing PET, PEESE, PET-PEESE, and Endogenous Kink (EK) methods
  * `waive()`: Robust extension with downweighting of spurious precision and outliers
* Features:
  * Instrumental variable approach using inverse sample sizes
  * Multiple weighting schemes (no weights, inverse-variance, MAIVE-adjusted, WAIVE)
  * Study-level correlation handling (fixed effects, clustering, or both)
  * Robust standard errors (CR0, CR1, CR2, wild bootstrap)
  * Anderson-Rubin confidence intervals for weak instruments
  * First-stage specification options (levels or log transformation)
  * Publication bias testing based on instrumented FAT
* Comprehensive test suite with 9 test files
* Documentation with examples and usage guidelines

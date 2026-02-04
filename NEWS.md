# MAIVE 0.2.4

*Released: 2026-02-04*

## Bug Fixes

* Easymeta.org links

---


# MAIVE 0.2.3

*Released: 2026-02-04*

## Internal

* Add warnings for weak instruments

---


# MAIVE 0.2.2

*Released: 2026-01-07*

## Internal

* Update the ar calculation to build the weighted residual correctly

---


# MAIVE 0.2.1

*Released: 2026-01-07*

## New Features

* Add the option to set an RNG seed at the highest function level

---


# MAIVE 0.2.0

*Released: 2026-01-07*

## New Features

* Add the ottawa conference slides link to strategic locations around the package
* Add an explicit citation file to the inst folder, update CLAUDE.md with instructions
* Add the funnel plot vignette; add vignette preview


## Bug Fixes

* Failing tests
* A couple of failing tests
* Makefile targets, update docs
* Add a missing funnel plot topic to pkgdown
* Funnel plot documentation


## Documentation

* Update the introduction vignette to include info on the available column mapping


## Internal

* Re-generate R docs
* Allow positional arguments in the main functions
* Add empty rows removal
* Add a validation module, update the main functions to utilize it
* Regenerate r docs
* Update references to WAIVE to highlight its more aggressive correction for phacking
* Update outdated paper references (2024 -> 2025, nature communications)
* Update the docs to feature links to easymeta.org
* Clean up more unused special characters
* Get rid of funnel_plot docstring special characters

---


# MAIVE 0.1.12

*Released: 2026-01-07*

## Bug Fixes

* Failing tests with a gt operator


## Internal

* Add an explicit function for generating a funnel plot

---


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

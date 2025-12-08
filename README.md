# MAIVE: Meta-Analysis Instrumental Variable Estimator

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/MAIVE)](https://CRAN.R-project.org/package=MAIVE)
[![R-CMD-check](https://github.com/meta-analysis-es/maive/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/meta-analysis-es/maive/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/meta-analysis-es/maive/branch/main/graph/badge.svg)](https://app.codecov.io/gh/meta-analysis-es/maive?branch=main)
<!-- badges: end -->

**Spurious Precision in Meta-Analysis of Observational Research**  
by Zuzana Irsova, Pedro R. D. Bom, Tomas Havranek, and Heiko Rachinger

**Project Website**: <https://meta-analysis.cz/maive/>

---

## Overview

MAIVE addresses a fundamental problem in meta-analysis of observational research: **spurious precision**.

Traditional meta-analysis assigns more weight to studies with lower standard errors, assuming higher precision. However, in observational research, precision can be manipulated through p-hacking and other questionable research practices, invalidating:

- Inverse-variance weighting schemes
- Traditional bias-correction methods (funnel plots, trim-and-fill)
- Selection models for publication bias

MAIVE implements an **instrumental variable approach** to limit bias caused by spurious precision in meta-analysis.

## Installation

### From CRAN (coming soon)

```r
install.packages("MAIVE")
```

### Development version

```r
install.packages("devtools")
devtools::install_github("meta-analysis-es/maive")
```

### Load package

```r
library(MAIVE)
```

## Quick Start

```r
# Prepare your data
data <- data.frame(
  bs = c(...),        # Effect sizes
  sebs = c(...),      # Standard errors
  Ns = c(...),        # Sample sizes
  study_id = c(...)   # Study IDs (optional)
)

# Run MAIVE with defaults (PET-PEESE, instrumented SEs, no weights)
result <- maive(
  dat = data,
  method = 3,      # PET-PEESE (default)
  weight = 0,      # No weights (default)
  instrument = 1,  # Instrument SEs (default)
  studylevel = 2,  # Cluster-robust (default)
  SE = 3,          # Wild bootstrap (default)
  AR = 1           # Anderson-Rubin CI (default)
)

# View results
print(result$Estimate)    # MAIVE estimate
print(result$SE)          # Standard error
print(result$Hausman)     # Hausman test
print(result$`F-test`)    # First-stage F-test
```

## Data Structure

The `maive()` function expects a data frame with:

| Column | Label | Description |
|--------|-------|-------------|
| 1 | `bs` | Primary estimates (effect sizes) |
| 2 | `sebs` | Standard errors (must be > 0) |
| 3 | `Ns` | Sample sizes (must be > 0) |
| 4 | `study_id` | Study identification (optional, for clustering/fixed effects) |

## Key Features

### Methods

- **PET** (Precision-Effect Test)
- **PEESE** (Precision-Effect Estimate with Standard Error)
- **PET-PEESE** (Conditional method, default)
- **EK** (Endogenous Kink)

### Weighting Schemes

- No weights (recommended when spurious precision is a concern)
- Inverse-variance weights
- MAIVE-adjusted weights (using instrumented SEs)
- **WAIVE** weights (robust downweighting via `waive()` function)

### Robust Inference

- Study-level correlation (fixed effects, clustering, or both)
- Multiple SE estimators (CR0, CR1, CR2, wild bootstrap)
- Anderson-Rubin confidence intervals for weak instruments
- First-stage specification options (levels or log transformation)

### Output

The function returns:

- MAIVE point estimate and standard error
- Standard (non-IV) estimate for comparison
- Hausman-type test statistic
- First-stage F-test of instrument strength
- Anderson-Rubin confidence interval
- Publication bias test p-value
- Instrumented standard errors

## Documentation

- **Getting Started Guide**: See `vignette("introduction")`
- **Function Reference**: `?maive` and `?waive`
- **Development Workflow**: `.github/DEVELOPMENT-WORKFLOW.md`
- **CRAN Submission**: `.github/CRAN-SUBMISSION.md`
- **Project Page**: <https://meta-analysis.cz/maive/>

## Example

```r
# Create example data
set.seed(123)
data <- data.frame(
  bs = rnorm(50, mean = 0.3, sd = 0.2),
  sebs = runif(50, min = 0.05, max = 0.3),
  Ns = sample(100:1000, 50, replace = TRUE),
  study_id = rep(1:10, each = 5)
)

# Run MAIVE
result <- maive(data, method = 3, weight = 0, instrument = 1, 
                studylevel = 2, SE = 3, AR = 1)

# Compare with standard estimate
cat("MAIVE Estimate:", result$Estimate, "\n")
cat("Standard Estimate:", result$StdEstimate, "\n")
cat("Hausman Test:", result$Hausman, "\n")

# Use WAIVE for robust estimation with outlier downweighting
result_waive <- waive(data, method = 3, instrument = 1, 
                      studylevel = 2, SE = 3, AR = 1)
cat("WAIVE Estimate:", result_waive$Estimate, "\n")
```

## Citation

If you use MAIVE in your research, please cite:

> Irsova, Z., Bom, P. R. D., Havranek, T., & Rachinger, H. (2024).
> Spurious Precision in Meta-Analysis of Observational Research.
> Available at: <https://meta-analysis.cz/maive/>

## References

Keane, M., & Neal, T. (2023). Instrument strength in IV estimation and inference: A guide to theory and practice. *Journal of Econometrics*, 235(2), 1625-1653. <https://doi.org/10.1016/j.jeconom.2022.12.009>

Tipton, E. (2015). Small sample adjustments for robust variance estimation with cluster-correlated data. *Psychological Methods*, 20(3), 375–389. <https://doi.org/10.1037/met0000019>

## Contributing

We welcome contributions! Please see our [GitHub repository](https://github.com/meta-analysis-es/maive) for:

- Bug reports and feature requests (use [Issues](https://github.com/meta-analysis-es/maive/issues))
- Code contributions (submit Pull Requests)
- Questions and discussions

## License

MIT License - see [LICENSE](LICENSE) file for details.

## Authors

- **Zuzana Irsova** - Charles University, Prague
- **Pedro R. D. Bom** - University of Deusto, Bilbao  
- **Tomas Havranek** - Charles University, Prague; CEPR, London; METRICS, Stanford
- **Petr Cala** (Maintainer) - Charles University, Prague
- **Heiko Rachinger** - University of the Balearic Islands, Palma

---

**Questions?** Contact the maintainer or visit our [project website](https://meta-analysis.cz/maive/).

# Test Suite

Run all package tests from the repository root:

```bash
Rscript -e 'pkgload::load_all(); testthat::test_dir("tests/testthat")'
```

This loads the MAIVE package in-place and executes the `testthat` specs. Ensure package dependencies are installed first (e.g. via `Rscript -e "devtools::install_dev_deps()"`).

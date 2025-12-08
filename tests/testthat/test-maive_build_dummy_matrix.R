test_that("maive_build_dummy_matrix produces one-hot encoded matrix", {
  values <- c("study_a", "study_b", "study_a", "study_c")
  result <- MAIVE:::maive_build_dummy_matrix(values)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), length(values))
  expect_equal(colnames(result), c("studyid.study_a", "studyid.study_b", "studyid.study_c"))
  expect_equal(result[, "studyid.study_a"], c(1, 0, 1, 0))
  expect_equal(result[, "studyid.study_b"], c(0, 1, 0, 0))
  expect_equal(result[, "studyid.study_c"], c(0, 0, 0, 1))
  expect_null(attr(result, "assign"))
  expect_null(attr(result, "contrasts"))
})

test_that("maive_build_dummy_matrix respects factor levels", {
  factor_values <- factor(c("x", "y", "x"), levels = c("y", "x"))
  result <- MAIVE:::maive_build_dummy_matrix(factor_values)

  expect_equal(colnames(result), c("studyid.y", "studyid.x"))
  expect_equal(result[, "studyid.y"], c(0, 1, 0))
  expect_equal(result[, "studyid.x"], c(1, 0, 1))
})

test_that("maive_build_dummy_matrix matches varhandle::to.dummy when available", {
  skip_if_not_installed("varhandle")

  library(varhandle)

  values <- c("A", "B", "A", "C", "B")
  expected <- varhandle::to.dummy(data.frame(studyid = values), "studyid")
  actual <- MAIVE:::maive_build_dummy_matrix(values)

  expect_identical(actual, expected)
})

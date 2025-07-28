test_that("coef_sd.pk() proceeds without error", {
  test_pk <- readRDS(paste0(test_path(), "/example_fitted_data.rds"))
  expect_no_error(coef_sd(test_pk))
})

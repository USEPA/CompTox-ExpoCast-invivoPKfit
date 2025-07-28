test_that("get_fit.pk proceeds without error", {
  test_pk <- readRDS(paste0(test_path(), "/example_fitted_data.rds"))
  expect_no_error(get_fit(test_pk))
  expect_s3_class(get_fit(test_pk), "data.frame")
})

test_that("get_data_original.pk proceeds without error", {
  test_pk <- readRDS(paste0(test_path(), "/example_fitted_data.rds"))
  expect_type(get_data_original(test_pk), "list")
  expect_s3_class(get_data_original(test_pk), class = "data.frame")
  expect_error(get_data_original(NULL))
})

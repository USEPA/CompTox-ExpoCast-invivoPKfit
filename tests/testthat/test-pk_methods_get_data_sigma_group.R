test_that("get_data_sigma_group.pk proceeds without error", {
  test_pk <- readRDS(paste0(test_path(), "/example_fitted_data.rds"))
  expect_no_error(get_data_sigma_group(test_pk))
  expect_type(get_data_sigma_group(test_pk), "integer")
  expect_s3_class(get_data_sigma_group(test_pk), "factor")
})


test_that("get_error_group.pk proceeds without error", {
  test_pk <- readRDS(paste0(test_path(), "/example_fitted_data.rds"))
  expect_no_error(get_error_group(test_pk))
  expect_no_error(get_error_group(test_pk, as_character = TRUE))
  expect_type(get_error_group(test_pk), "list")
  expect_type(get_error_group(test_pk, as_character = TRUE), "character")
})

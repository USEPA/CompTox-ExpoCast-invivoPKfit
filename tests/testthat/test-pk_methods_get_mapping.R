test_that("get_mapping.pk proceeds without error", {
  test_pk <- readRDS(paste0(test_path(), "/example_fitted_data.rds"))
  expect_no_error(get_mapping(test_pk))
  expect_s3_class(get_mapping(test_pk), "uneval")
})

test_that("get_data_sigma_group.pk proceeds without error", {
  test_pk <- readRDS(paste0(test_path(), "/example_fitted_data.rds"))
  expect_no_error(get_status(test_pk))
  expect_equal(get_status(test_pk), 5L)
  expect_warning(test_pk <- do_preprocess(test_pk))
  expect_no_error(get_status(test_pk))
  expect_equal(get_status(test_pk), 2L)
})

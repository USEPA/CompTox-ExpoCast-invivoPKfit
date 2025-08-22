test_that("get_data_info.pk proceeds without error", {
  test_pk <- readRDS(paste0(test_path(), "/example_fitted_data.rds"))
  expect_no_error(get_data_info(test_pk))
  expect_type(get_data_info(test_pk), "list")
  expect_true(all(names(get_data_info(test_pk)) %in% c("data_summary", "data_flags", "dose_norm_check", "nca")))
})

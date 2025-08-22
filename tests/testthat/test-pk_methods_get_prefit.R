test_that("get_prefit.pk proceeds without error", {
  test_pk <- readRDS(paste0(test_path(), "/example_fitted_data.rds"))
  expect_no_error(get_prefit(test_pk))
  expect_type(get_prefit(test_pk), "list")
  expect_true(all(c("sigma_DF", "par_DF", "fit_check") %in% names(get_prefit(test_pk))))
})

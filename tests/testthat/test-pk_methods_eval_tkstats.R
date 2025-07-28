test_that("logLike.pk proceeds without error", {
  test_pk <- readRDS(paste0(test_path(), "/example_fitted_data.rds"))
  library(httk)
  expect_no_error(eval_tkstats(test_pk))
})

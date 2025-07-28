test_that("AIC.pk proceeds without error", {
  test_pk <- readRDS(paste0(test_path(), "/example_fitted_data.rds"))
  expect_no_error(AIC(test_pk))
})

test_that("AIC.pk proceeds without error with newdata argument", {
  test_pk <- readRDS(paste0(test_path(), "/example_fitted_data.rds"))
  test_data <- get_data.pk(test_pk) |>
    dplyr::group_by(DATA_GROUP_ID) |>
    dplyr::slice_head(prop = 0.9) |>
    dplyr::ungroup()
  expect_no_error(AIC(test_pk, newdata = test_data))
})

test_that("BIC.pk proceeds without error", {
  test_pk <- readRDS(paste0(test_path(), "/example_fitted_data.rds"))
  expect_no_error(BIC(test_pk))
})

test_that("BIC.pk proceeds without error with newdata argument", {
  test_pk <- readRDS(paste0(test_path(), "/example_fitted_data.rds"))
  test_data <- get_data.pk(test_pk) |>
    dplyr::group_by(DATA_GROUP_ID) |>
    dplyr::slice_head(prop = 0.9) |>
    dplyr::ungroup()
  expect_no_error(BIC(test_pk, newdata = test_data))
})


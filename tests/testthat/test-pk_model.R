test_that("pk_model constructor works", {
  expect_no_error(pk_model(
    name = "model_1comp",
    params = c(
      "kelim", "Vdist", "Fgutabs",
      "kgutabs", "Fgutabs_Vdist", "Rblood2plasma"
    ),
    conc_fun = "cp_1comp",
    auc_fun = "auc_1comp",
    params_fun = "get_params_1comp",
    tkstats_fun = "tkstats_1comp",
    conc_fun_args = NULL,
    auc_fun_args = NULL,
    params_fun_args = NULL,
    tkstats_fun_args = NULL
  ))
})

test_that("Checking pk_model object class works", {
  expect_true(is.pk_model(model_1comp))
  expect_false(is.pk_model("model_1comp"))
})

test_that("Setting multiple parameters to optimize works!", {
  expect_no_error(
    test_model <- set_params_optimize(model_1comp, params = c("kelim", "kgutabs"))
  )
  expect_identical(test_model$params_fun_args$pars_to_optimize, c("kelim", "kgutabs"))
})

test_that("Setting parameter starts in `model$param_fun_args` works!", {
  expect_no_error(
    test_model <- set_params_starts(model_2comp, starts = list(k12 = 0.45))
  )
  expect_identical(test_model$params_fun_args$param_starts, list(k12 = 0.45))
})


test_that("Empty parameters argument errors when trying to set to optimize", {
  expect_error(set_params_optimize(model_1comp, params = c("k12", "k21")))
})

test_that("Non-model parameters are discarded with a warning but valid parameters are kept", {
  expect_warning(
    test_model <- set_params_optimize(
      model_1comp,
      params = c("kelim", "k12", "k21")
    )
  )
  expect_identical(test_model$params_fun_args$pars_to_optimize, "kelim")
})

test_that("Setting using defaults and using parameter groups works!", {
  expect_no_error(
    test_model <- pk_model(
      name = "model_1comp",
      params = c(
        "kelim", "Vdist", "Fgutabs",
        "kgutabs", "Fgutabs_Vdist", "Rblood2plasma"
      ),
      conc_fun = "cp_1comp",
      auc_fun = "auc_1comp",
      params_fun = "get_params_1comp",
      tkstats_fun = "tkstats_1comp",
      conc_fun_args = NULL,
      auc_fun_args = NULL,
      params_fun_args = NULL,
      tkstats_fun_args = NULL,
      param_groups = list(default = c("kelim", "Vdist", "Fgutabs", "kgutabs", "Rblood2plasma"))
    )
  )
  expect_no_error(test_model <- set_params_optimize(test_model))
  expect_true(
    all(test_model$param_fun_args$pars_to_optimize %in%
      c("kelim", "Vdist", "Fgutabs", "kgutabs", "Rblood2plasma"))
  )
})


test_that("Toggling clearance parameter in `pk_models` works!", {
  expect_no_error(test_model <- toggle_clearance_mode(model_httk_gas_pbtk))
  expect_false(test_model$conc_fun_args$restrictive)
  expect_message(test_model <- toggle_clearance_mode(model_httk_gas_pbtk))
  expect_no_error(toggle_clearance_mode(model_1comp))
})

test_that(
  "Checking model parameters with both parameter groups and extra parameter gives warning",
  {
    expect_no_error(test_model <- pk_model(
      name = "model_1comp",
      params = c(
        "kelim", "Vdist", "Fgutabs",
        "kgutabs", "Fgutabs_Vdist", "Rblood2plasma"
      ),
      conc_fun = "cp_1comp",
      auc_fun = "auc_1comp",
      params_fun = "get_params_1comp",
      tkstats_fun = "tkstats_1comp",
      conc_fun_args = NULL,
      auc_fun_args = NULL,
      params_fun_args = NULL,
      tkstats_fun_args = NULL,
      param_groups = list(default = c("kelim", "Vdist", "Fgutabs"))
    ))
    expect_warning(set_params_optimize(test_model, params = c("kelim", "default")))
  }
)



test_that(
  "Checking that bespoke models with constant set starts can be fit",
  {
    test_bespoke <- function() {
      test_model_1comp <- set_params_optimize(
        model = model_1comp,
        params = c("kelim", "kgutabs")
      ) |>
        set_params_starts(starts = list(Vdist = 1))
      test_model_1comp <- adjust_model_name(test_model_1comp)

      test_pk <- pk(
        data = subset(
          cvt,
          analyzed_chem_dtxsid %in% c("DTXSID0020442", "DTXSID0021961") & species %in% "rat"
        )
      ) + stat_model(c("model_flat", "test_model_1comp"))

      expect_no_error(
        test_pk <- do_fit.pk(test_pk)
      )

      return(test_pk)
    }

    expect_no_error(test_fit <- test_bespoke())
    expect_no_error(
      test_fit <- subset(test_fit$fit, model == "test_model_1comp" & param_name == "Vdist")
    )

    expect_equal(unique(unlist(test_fit$estimate)), 1.0)

    expect_equal(
      unique(unlist(test_fit$start)),
      unique(unlist(test_fit$estimate))
    )
  }
)

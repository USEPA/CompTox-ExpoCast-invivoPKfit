test_that(
  "fitting possible using 'auto' time units proceed without error",
  code = {
    test_pk <- pk(
      data = subset(cvt, analyzed_chem_dtxsid %in% c("DTXSID5048265", "DTXSID3031860"))
    ) + scale_time(new_units = "auto") +
      settings_preprocess(suppress.messages = TRUE)
    expect_no_error(test_pk <- do_fit(test_pk))
  }
)





test_that(
  "fitting possible using concentration scaling proceed without error",
  {
    test_pk <- pk(
      data = subset(
        cvt,
        analyzed_chem_dtxsid %in% "DTXSID1021116" & species %in% c("rat", "human")
      )
    ) + scale_conc(dose_norm = TRUE, log10_trans = TRUE) +
      settings_preprocess(suppress.messages = TRUE)
    expect_no_error(test_pk <- do_fit(test_pk))
  }
)

test_that(
  "fitting for all models (1-compartment, 2-compartment, flat, and gas_pbtk) is possible",
  {
    test_pk <- pk(
      data = subset(
        cvt,
        analyzed_chem_dtxsid %in% c("DTXSID0020232", "DTXSID2023309") & species %in% "rat"
      )
    ) +
      stat_model(model = c("model_flat", "model_1comp", "model_2comp", "model_httk_gas_pbtk")) +
      settings_preprocess(suppress.messages = TRUE)
    expect_no_error(test_pk <- do_fit(test_pk))
  }
)




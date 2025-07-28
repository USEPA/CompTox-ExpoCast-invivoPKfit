test_that(
  "fitting possible using 'auto' time units proceed without error",
  code = {
    test_pk <- pk(
      data = subset(cvt, analyzed_chem_dtxsid %in% "DTXSID5048265")
    ) + scale_time(new_units = "auto")
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
    ) + scale_conc(dose_norm = TRUE, log10_trans = TRUE)
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
    )
    expect_no_error(test_pk <- do_fit(test_pk))
  }
)




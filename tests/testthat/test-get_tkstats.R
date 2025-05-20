test_that("get_tkstats works without error when expected",
          {
        test_pk <- pk(data = subset(
          cvt,
          analyzed_chem_dtxsid %in% "DTXSID1021116" & species %in% "rat"
          ))
        test_pk <- test_pk + stat_model()
        test_pk <- do_fit(test_pk)
        expect_no_error(get_tkstats(test_pk))
          })

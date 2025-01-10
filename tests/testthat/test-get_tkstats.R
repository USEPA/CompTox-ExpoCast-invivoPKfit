test_that("get_tkstats works without error when expected",
          {
        my_pk <- pk(data = subset(
          cvt,
          analyte_dtxsid %in% "DTXSID1021116" & species %in% "rat"
          ))
        my_pk <- my_pk + stat_model()
        my_pk <- do_fit(my_pk)
        expect_no_error(get_tkstats(my_pk))
          })

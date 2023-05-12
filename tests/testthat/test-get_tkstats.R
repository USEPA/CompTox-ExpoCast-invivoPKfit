test_that("get_tkstats works without error when expected",
          {
        my_pk <- pk(data = subset(
          cvt,
          chemicals_analyzed.dsstox_substance_id %in% "DTXSID1021116" &
            subjects.species %in% "rat"
          ))
        my_pk <- my_pk + stat_model()
        my_pk <- fit(my_pk)
        expect_no_error(get_tkstats(my_pk))
          })

test_that("data preprocessing works without errors",
          {
            my_pk <- pk(
              data = subset(cvt,
                            analyte_dtxsid %in% "DTXSID3061635"
              )
            )
            expect_no_error(do_preprocess(my_pk))
          }
)

test_that("data preprocessing handles NULL data as expected",
          {
            my_pk <- pk(
              data = NULL
              )
            expect_message(do_preprocess(my_pk),
                            regexp = "do_preprocess.pk(): Original data is NULL",
                            fixed = TRUE)
          }
)

test_that("data preprocessing adds an element 'data'",
          {
            my_pk <- pk(
              data = subset(cvt,
                            analyte_dtxsid %in% "DTXSID3061635"
              ))
            my_pk <- do_preprocess(my_pk)
            expect_true("data" %in% names(my_pk))
          }
)

test_that("preprocessed data has all of the required harmonized variable names",
          {
            my_pk <- pk(
              data = subset(cvt,
                            analyte_dtxsid %in% "DTXSID3061635"
              )
            )
            my_pk <- do_preprocess(my_pk)
            required_vars <- c("Chemical",
                               "Species",
                               "Route",
                               "Media",
                               "Time",
                               "Time.Units",
                               "Dose",
                               "Dose.Units",
                               "Value",
                               "Value.Units",
                               "LOQ",
                               "Value_SD",
                               "Conc",
                               "Conc.Units",
                               "Detect",
                               "Conc_SD",
                               "N_Subjects",
                               "Time_trans",
                               "Time_trans.Units",
                               "Conc_trans",
                               "Conc_trans.Units",
                               "Conc_SD_trans")
            expect_true(
              all(
                required_vars %in% names(my_pk$data)
              )
            )
          }
)

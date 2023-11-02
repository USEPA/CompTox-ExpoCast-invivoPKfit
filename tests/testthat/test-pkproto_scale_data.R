test_that("scale_conc() returns an object of class pkproto",
          {
            expect_true(is.pkproto(scale_conc()))
          }
)

test_that("scale_conc() returns an object of class pk_scales",
          {
            expect_true(is.pk_scales(scale_conc()))
          }
)

test_that("scale_conc() returns an object with element 'name' = 'conc'",
          {
            obj <- scale_conc()
            expect_equal(obj$name, "conc")
          }
)

test_that("scale_conc() returns an object with the expected names in its element 'value'",
          {
            obj <- scale_conc()
            expect_true(
              all(
                c("dose_norm",
                  "log10_trans",
                  "expr") %in%
                  names(obj$value)
              )
            )
          }
)


test_that("scale_time() returns an object of class pkproto",
          {
            expect_true(is.pkproto(scale_time()))
          }
)

test_that("scale_time() returns an object of class pk_scales",
          {
            expect_true(is.pk_scales(scale_time()))
          }
)

test_that("scale_time() returns an object with element 'name' = 'time'",
          {
            obj <- scale_time()
            expect_equal(obj$name, "time")
          }
)

test_that("scale_time() returns an object with the expected names in its element 'value'",
          {
            obj <- scale_time()
            expect_true(
              all(
                c("new_units") %in%
                  names(obj$value)
              )
            )
          }
)

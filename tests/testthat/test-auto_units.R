test_that(
  "time points over 2 days get converted to hours",
  {
    test_time <- seq(0,2,by=(1/24))
    expect_equal(auto_units(y = test_time,
                            from = "days"),
                 "hours")
  }
)

test_that(
  "time points over 24 hours stay in hours",
  {
    test_time <- seq(0,24)
    expect_equal(auto_units(y = test_time,
                            from = "hours"),
                 "hours")
  }
)

test_that(
  "time points in hours over a month get converted to days",
  {
    expect_equal(
      auto_units(y = convert_time(seq(0,30),
                                  from = "days",
                                  to = "hours"),
                 from = "hours"),
      "days"
    )
  }
)

test_that(
  "time points in hours over a year get converted to months",
  {
    expect_equal(
      auto_units(y = convert_time(seq(0,365),
                                  from = "days",
                                  to = "hours"),
                 from = "hours"),
      "months"
    )
  }
)

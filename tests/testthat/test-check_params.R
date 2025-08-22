test_that("check_params_flat errors when missing relevant parameters", {
  expect_identical(
    check_params_flat(params = list(Rblood2plasma = 0.86), route = "iv", medium = "blood"),
    "For flat/NULL IV model, missing parameters: Vdist"
  )
  expect_identical(
    check_params_flat(params = list(Rblood2plasma = 0.86), route = "oral", medium = "blood"),
    "For flat/NULL oral model, missing parameters: Fgutabs_Vdist"
  )
})


test_that("check_params_1comp errors when missing relevant parameters", {
  expect_identical(
    check_params_1comp(params = list(Rblood2plasma = 0.86,
                                    kelim = 0.2), route = "iv", medium = "blood"),
    "For 1-compartment IV model, missing parameters: Vdist"
  )
  expect_identical(
    check_params_1comp(params = list(Rblood2plasma = 0.86,
                                    kelim = 0.2, Fgutabs_Vdist = 2.0), route = "oral", medium = "blood"),
    "For 1-compartment oral model, missing parameters: kgutabs"
  )
})

test_that("check_params_2comp errors when missing relevant parameters", {
  expect_identical(
    check_params_2comp(params = list(Rblood2plasma = 0.86,
                                    kelim = 0.2), route = "iv", medium = "blood"),
    "For 2-compartment IV model, missing parameters: V1, k12, and k21"
  )
  expect_identical(
    check_params_2comp(params = list(Rblood2plasma = 0.86,
                                     kelim = 0.2, k12 = 2.0, k21 = 0.5), route = "oral", medium = "blood"),
    "For 2-compartment oral model, missing parameters: Fgutabs_V1 and kgutabs"
  )
})

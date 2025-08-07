devtools::load_all()

#' 1-compartment model object.
model_1comp <- pk_model(
  name = "model_1comp",
  params = c(
    "kelim",
    "Vdist",
    "Fgutabs",
    "kgutabs",
    "Fgutabs_Vdist",
    "Rblood2plasma"
  ),
  conc_fun = "cp_1comp",
  auc_fun = "auc_1comp",
  params_fun = "get_params_1comp",
  tkstats_fun = "tkstats_1comp",
  conc_fun_args = NULL,
  auc_fun_args = NULL,
  params_fun_args = NULL,
  param_groups = list(
    default = c(
      "kelim",
      "Vdist",
      "Fgutabs",
      "kgutabs",
      "Fgutabs_Vdist",
      "Rblood2plasma"
    ),
    rates = c("kelim", "kgutabs")
  ),
  tkstats_fun_args = NULL
) |> set_params_optimize()


model_2comp <- pk_model(
  name = "model_2comp",
  params = c(
    "kelim",
    "k12",
    "k21",
    "V1",
    "Fgutabs",
    "kgutabs",
    "Fgutabs_V1",
    "Rblood2plasma"
  ),
  conc_fun = "cp_2comp",
  auc_fun = "auc_2comp",
  params_fun = "get_params_2comp",
  tkstats_fun = "tkstats_2comp",
  conc_fun_args = NULL,
  auc_fun_args = NULL,
  params_fun_args = NULL,
  param_groups = list(
    default = c(
      "kelim",
      "k12",
      "k21",
      "V1",
      "Fgutabs",
      "kgutabs",
      "Fgutabs_V1",
      "Rblood2plasma"
    ),
    rates = c("kelim", "k12", "k21", "kgutabs")
  ),
  tkstats_fun_args = NULL
) |> set_params_optimize()

model_flat <- pk_model(
  name = "model_flat",
  params = c(
    "Vdist",
    "Fgutabs",
    "Fgutabs_Vdist",
    "Rblood2plasma"
  ),
  conc_fun = "cp_flat",
  auc_fun = "auc_flat",
  params_fun = "get_params_flat",
  tkstats_fun = "tkstats_flat",
  conc_fun_args = NULL,
  auc_fun_args = NULL,
  params_fun_args = NULL,
  param_groups = list(
    default = c(
      "Vdist",
      "Fgutabs",
      "Fgutabs_Vdist",
      "Rblood2plasma"
    )
  ),
  tkstats_fun_args = NULL
) |> set_params_optimize()

model_httk_gas_pbtk <- pk_model(
  name = "model_httk_gas_pbtk",
  params = c(
    "BW",
    "Caco2.Pab",
    "Caco2.Pab.dist",
    "Clint",
    "Clint.dist",
    "Clmetabolismc",
    "Funbound.plasma",
    "Funbound.plasma.dist",
    "Funbound.plasma.adjustment",
    "Fabsgut",
    "Fhep.assay.correction",
    "hematocrit",
    "Kgut2pu",
    "Krbc2pu",
    "kgutabs",
    "Kkidney2pu",
    "Klung2pu",
    "km",
    "Kmuc2air",
    "Kliver2pu",
    "Krest2pu",
    "Kblood2air",
    "kUrtc",
    "liver.density",
    "logHenry",
    "million.cells.per.gliver",
    "MW",
    "Pow",
    "pKa_Donor",
    "pKa_Accept",
    "MA",
    "Qcardiacc",
    "Qgfrc",
    "Qgutf",
    "Qliverf",
    "Qalvc",
    "Qkidneyf",
    "Qlungf",
    "Rblood2plasma",
    "Vgutc",
    "Vliverc",
    "Vartc",
    "Vkidneyc",
    "Vlungc",
    "vmax",
    "Vmucc",
    "Vvenc",
    "Vrestc",
    "KFsummary",
    "Fprotein.plasma",
    "fabs.oral",
    "Qgut",
    "Qintesttransport"
  ),
  conc_fun = "cp_httk_gas_pbtk",
  auc_fun = "auc_httk_gas_pbtk",
  params_fun = "get_params_httk_gas_pbtk",
  tkstats_fun = "tkstats_httk_gas_pbtk",
  conc_fun_args = alist(
    restrictive = TRUE,
    this_chem = Chemical,
    this_species = Species
  ),
  auc_fun_args = NULL,
  params_fun_args = list(restrictive = TRUE),
  param_groups = list(
    default = c("Clint", "Funbound.plasma"),
    hepatic_clearance = c("Clint", "Funbound.plasma"),
    partition_coefficients = c(
      "Kgut2pu",
      "Krbc2pu",
      "Kkidney2pu",
      "Klung2pu",
      "Kmuc2air",
      "Kliver2pu",
      "Krest2pu",
      "Kblood2air"
    )
  ),
  tkstats_fun_args = alist(
    restrictive = TRUE,
    this_chem = Chemical,
    this_species = Species
  )
) |> set_params_optimize()


usethis::use_data(model_1comp,
  model_2comp,
  model_flat,
  model_httk_gas_pbtk,
  internal = FALSE,
  overwrite = TRUE
)

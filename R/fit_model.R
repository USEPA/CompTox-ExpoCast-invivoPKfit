fit_model <- function(data,
                      model,
                      model_fun,
                      ratio_data_to_dose,

                      compound_col = "compound",
                      compound_default = NULL,

                      cas_col = "cas",
                      cas_default = NULL,

                      reference_col = "reference",
                      reference_default = NULL,

                      species_col = "species",
                      species_default = NULL,

                      species_weight_col = "species_weight",
                      species_weight_default = NULL,
                      species_weight_units_col = "species_weight_units",
                      species_weight_units_default = NULL,

                      dose_col = "dose",
                      dose_default = NULL,

                      time_col = "time",
                      time_default = NULL,
                      time_units_col = "time_units",
                      time_units_default = NULL,

                      media_col = "media",
                      media_default = NULL,
                      media_units_col = "media_units",
                      media_units_default = NULL,

                      conc_col = "conc",
                      conc_default = NULL,
                      conc_units_col = "conc_units",
                      conc_units_default = NULL,

                      route_col = "route",
                      route_default = NULL,

                      loq_col = "loq",
                      loq_default = NULL,

                      loq_units_col = "loq_units",
                      loq_units_default = NULL) {

  data <- rename_columns(data,

                         compound_col,
                         compound_default,

                         cas_col,
                         cas_default,


                         reference_col,
                         reference_default,

                         species_col,
                         species_default,

                         species_weight_col,
                         species_weight_default,
                         species_weight_units_col,
                         species_weight_units_default,

                         dose_col,
                         dose_default,

                         time_col,
                         time_default,
                         time_units_col,
                         time_units_default,

                         media_col,
                         media_default,
                         media_units_col,
                         media_units_default,

                         conc_col,
                         conc_default,
                         conc_units_col,
                         conc_units_default,

                         route_col,
                         route_default,

                         loq_col,
                         loq_default,

                         loq_units_col,
                         loq_units_default)

  ### convert to data.table
  data <- data.table::as.data.table(data)

  data <- clean_data(data, ratio_data_to_dose)

  if (model == "noncompartment") {

    pk_fit_table <- fit_noncomp(data)

  } else {

    params_by_cas_spec <- data[, unique(.SD[, .(compound)]), by = .(cas, species)]

  }

  if (model == "1compartment") {

    ### Wambaugh et al. (2018) medians
    params_by_cas_spec[, kelim := 0.25]
    params_by_cas_spec[, vdist := 5.56]
    params_by_cas_spec[, fgutabs := 1.0]
    params_by_cas_spec[, kgutabs := 2.19]

  } else if (model == "2compartment") {

    if (any(params_by_cas_spec$cas %in% httk::get_cheminfo())) {

      params_by_cas_spec[cas %in% httk::get_cheminfo(),
                         c("kelim",
                           "vdist",
                           "fgutabs",
                           "kgutabs") := httk::parameterize_1comp(chem.cas = cas,
                                                                  default.to.human = TRUE,
                                                                  species = species)[c("kelim",
                                                                                       "vdist",
                                                                                       "fgutabs",
                                                                                       "kgutabs")],
                         by = c("cas", "species")]
    }

  }


  return(pk_fit_table)
}

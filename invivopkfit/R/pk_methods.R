#'Create a new `pk` object
#'
#'
#'[pk()] initializes a new `pk` object. It can be used to declare the input data
#'frame for a PK analysis.
#'
#'
#'[pk()] is used to construct the initial `pk` object for analysis, i.e., to
#'provide the concentration-dose-time data to be fitted. It is almost always
#'followed by `+` to add steps to the workflow.
#'
#'For example, you could use `pk(my_data) + stat_model(model = '1comp')` to set
#'up and fit a 1-compartment model.
#'
#'# Mappings
#'
#'Your input data can have any variable names you like. However, internally,
#'`invivopkfit` needs to use a set of "standard", harmonized variable names
#'(e.g., it refers to the variable containing measured tissue concentrations as
#'`Conc`). The `mapping` argument sets up a mapping between the variable names
#'in `data`, and the internal harmonized variable names used by `invivopkfit`.
#'
#'The simplest, most flexible way to set up this mapping is by (ab)using a call
#'to [ggplot2::aes()]. Usually, you would use [ggplot2::aes()] to set up
#'mappings to aesthetic variables for plotting: `x`, `y`, `color`, `size`,
#'`shape`, and so on. Here, instead, we are setting up mappings to `invivopkfit`
#'variables used in model fitting. These "`invivopkfit` aesthetic" variables are
#'as follows:
#'
#'
#' - `DTXSID`: A `character` variable containing the chemical's DSSTox ID. All rows of `data` should have the same value for `DTXSID`.
#' - `Species`: A `character` variable containing the name of the species for which the data was measured.  All rows of `data` should have the same value for `Species`.
#' - `Reference`: A `character` variable containing a unique identifier for the data reference (e.g., a single publication).
#' - `Subject`: A `character` variable containing a unique identifier for the subject associated with each observation (an individual animal or group of animals).
#' - `N_Subjects`: A `numeric` variable; an integer giving the number of individual animals represented by this observation. (Some data sets report tissue concentrations for individual animals, in which case `N_Subjects` will be 1; others report average tissue concentrations for groups of multiple animals, in which case `N_Subjects` will be greater than 1.)
#' - `Weight`: A `numeric` variable giving the subject's body weight.
#' - `Weight.Units`: A `character` variable giving the units of body weight.
#' - `Route`: A `character` variable denoting the route of administration. Either `po` (oral) or `iv` (intravenous). Other routes are not currently supported.
#' - `Dose`: A `numeric` variable giving the dose administered.
#' - `Dose.Units`: A `character` variable giving the units of the administered doses.
#' - `Time`: A `numeric` variable giving the time of tissue collection.
#' - `Time.Units`: A `numeric` variable giving the units of `Time`.
#' - `Media`: A `character` variable giving the tissue that was analyzed. Either `blood` or `plasma`. Other tissues are not currently supported.
#' - `Value`: A `numeric` variable giving the tissue concentration in units of mg/L. If `N_Subjects > 1`, `Value` is assumed to represent the mean tissue concentration for this group of subjects. If the tissue concentration was below the limit of quantification (LOQ), this value may be `NA_real_`.
#' - `Value_SD`: A `numeric` variable giving the standard deviation of the tissue concentration in units of mg/L, if available and relevant. If `N_Subjects > 1`, `Value_SD` is assumed to represent the standard deviation of tissue concentrations for this group of subjects. If `N_Subjects == 1`, then `Value_SD` may be `NA_real_`.
#' - `LOQ`: A `numeric` variable giving the limit of quantification applicable to this tissue concentration in units of mg/L, if available.
#' - `Value.Units`: A `character` variable giving the units of `Value`, `Value_SD`, and `LOQ`.
#'
#'As with usual calls to [ggplot2::aes()], you should provide the variable names
#'without quoting them. For example, use `ggplot2::aes(DTXSID = my_dtxsid)`. Do
#'*not* use `ggplot2::aes("DTXSID" = "my_dtxsid")`.
#'
#'Also, as with usual calls to [ggplot2::aes()], you may also specify that any
#'of the "`invivopkfit` aesthetic" variables should be mapped to a constant
#'value, rather than to a variable in `data`. For example, imagine that you
#'don't have a column in `data` that encodes the units of body weight, but you
#'know that all body weights are provided in units of kilograms. You could
#'specify `mapping = ggplot2::aes(DTXSID = my_dtxsid, Species = my_species,
#'Weight = my_weight, Weight.Units = "kg")` to map `Weight.Units` to a fixed
#'value of "kg".
#'
#'Finally, as with usual calls to [ggplot2::aes()], you may specify mappings as
#'expressions that use variable names in `data`. For example, if the
#'species-name variable in `data` sometimes says "rat", sometimes "Rat",
#'sometimes "RAT", you might want to harmonize the capitalization. You can do
#'that easily by specifying `mapping = ggplot2::aes(DTXSID = my_dtxsid, Species
#'= tolower(my_species)`.
#'
#'# Data
#'
#'`data` should contain data for only one `DTXSID` and one `Species`. It may
#'contain data for multiple `Route` and/or `Media` values, which can be fitted
#'simultaneously. It may contain data from multiple `Study` IDs (in the case
#'where a joint or pooled analysis is desired), or from only one `Study` ID (in
#'the case where a separate analysis is desired for one reference at a time).
#'
#'# Data settings
#'
#'The default value of argument `data_settings` is a named list with the
#'following elements:
#'- `ratio_conc_to_dose` Ratio between the mass units used to report the
#'concentration data and the mass units used to report the dose. Default 1. For
#'example, concentration reported in ug/L and dose reported in mg/kg/day would
#'require `ratio_conc_to_dose = 0.001`, because 1 ug/1 mg = 1e-6 g / 1e-3 g =
#'0.001.
#' - `routes_keep` List of routes to retain. Default `c("po", "iv")` to
#'retain only oral and IV administration data.
#'- `media_keep` List of media to retain. Default `c("blood", "plasma")` to
#'retain only concentrations in blood and plasma.
#' - `impute_loq` Logical: TRUE to impute values for missing LOQs; FALSE to
#'leave them alone.
#' - `impute_sd` Logical: TRUE to impute values for missing sample SDs for
#'multi-subject observations; FALSE to leave them alone
#' - `suppress.messages` Logical: Whether to suppress verbose messages.
#'Default FALSE, to be verbose.
#'
#'# Scales
#'
#'The optional argument `scales` is a named list with two elements:
#'
#' - `conc` is itself a list with two elements, `normalize` and `trans`. By default, both are `"identity"`. See [scale_conc()] for options.
#' - `time` is itself a list with one element, `trans`. By default, it is `"identity"`. See [scale_time()] for options.
#'
#'You can always modify `scales` later by using [scale_conc()] and/or
#'[scale_trans()], and it may actually be easier to do it that way. For example,
#'`pk(data = my_df, scales = list(conc = list(trans = log10)))` does exactly the
#'same thing as `pk(data = my_df) + scale_conc(trans = log10)`.
#'
#'# optimx settings
#'
#'Argument `optimx_settings` is a named list containing settings for the
#'optimizer, with the following elements:
#'
#'  - `method` The method to use, from those implemented in [optimx::optimx()]. Default `"bobyqa"`.
#' - `itnmax` The maximum number of iterations, as in [optimx::optimx()]. Default `1e6`.
#' - `hessian` Whether to compute the Hessian after optimizing, as in [optimx::optimx()]. Default `FALSE`.
#' - `control` A named list of control parameters for the optimizer. See [optimx::optimx()] for options and defaults.  Default here is `list(kkt = FALSE)`.
#'
#'
#'
#'@param data A `data.frame`. The default is an empty data frame.
#'@param mapping A mapping set up by (ab)using [ggplot2::aes()]. Call is of form
#'  `ggplot2::aes(new_variable = old_variable)` `new_variable` represents the
#'  harmonized variable name that will be used within `invivopkfit`;
#'  `old_variable` represents the variable name in `data`. If you want to
#'  provide a fixed/constant value for a `new_variable` rather than taking its
#'  value from a variable in `data`, simply supply that fixed/constant value in
#'  the `old_variable` position.
#'@param data_settings A named list of settings for the data preprocessing. See
#'  Details for options and defaults.
#'@param scales A named list of scales (normalization and/or transformations) to
#'  use for concentration and time data. See Details for options and defaults.
#'@param optimx_settings A named list of settings for the optimizer. See Details
#'  for options and defaults.
#'@return An object of class `pk`. The initial `pk` object is a list with
#'  elements `data_orig`, `data_settings`, `scales` and `optimx_settings`.
#'  `data_orig` is the original data set to be fitted, as supplied in the
#'  argument `data`. `data_settings` is a named list containing all the other
#'  input arguments: these provide settings that will be used when the data is
#'  pre-processed before fitting.
#'@author Caroline Ring
#'@export

pk <- function(data = NULL,
               mapping = ggplot2::aes(
                 Compound_Dosed = studies.test_substance_name_original,
                 DTXSID_Dosed = chemicals_dosed.dsstox_substance_id,
                 CAS_Dosed = chemicals_dosed.dsstox_casrn,
                 Compound_Analyzed = series.analyte_name_original,
                 DTXSID_Analyzed = chemicals_analyzed.dsstox_substance_id,
                 CAS_Analyzed = chemicals_analyzed.dsstox_casrn,
                 Reference = documents_reference.id,
                 Extraction = documents_extraction.id,
                 Species = subjects.species,
                 Weight =subjects.weight_kg,
                 Weight.Units = "kg",
                 Dose = studies.dose_level_normalized_corrected,
                 Dose.Units = "mg/kg",
                 Time = conc_time_values.time_hr,
                 Time.Units = "hours",
                 Media = series.conc_medium_normalized,
                 Value = conc_time_values.conc,
                 Value_SD = conc_time_values.conc_sd_normalized,
                 Value.Units = "mg/L",
                 Route = studies.administration_route_normalized,
                 LOQ = series.loq_normalized,
                 Subject = subjects.id,
                 N_Subjects = series.n_subjects_in_series,
                 Study_ID = studies.id,
                 Series_ID = series.id
               ),
               data_settings = list(
                 ratio_conc_to_dose = 1,
               calc_loq_factor = 0.45,
               routes_keep = c("po", "iv"),
               media_keep = c("blood", "plasma"),
               impute_loq = TRUE,
               impute_sd = TRUE,
               suppress.messages = FALSE
               ),
               scales = list(conc = list(normalize = "identity",
                                           trans = "identity"),
                               time = list(trans = "identity")),
               optimx_settings = list(
                 method = "bobyqa",
               itnmax = 1e6,
               hessian = FALSE,
               control = list(kkt = FALSE)
               )
               ){

  #fill in any non-specified data settings with their defaults
  data_settings_default <- list(
    ratio_conc_to_dose = 1,
    calc_loq_factor = 0.45,
    routes_keep = c("po", "iv"),
    media_keep = c("blood", "plasma"),
    impute_loq = TRUE,
    impute_sd = TRUE,
    suppress.messages = FALSE
  )

  missing_data_settings <- setdiff(names(data_settings_default),
                                   names(data_settings))
  data_settings[missing_data_settings] <- data_settings_default[missing_data_settings]

  #fill in any non-specified scales with their defaults
  scales_default <- list(conc = list(normalize = "identity",
                   trans = "identity"),
       time = list(trans = "identity"))

  missing_scales <- setdiff(names(scales_default),
                            names(scales))
  scales[missing_scales] <- scales_default[missing_scales]
#check for missing scales_conc
  missing_scales_conc <- setdiff(names(scales_default$conc),
                                 names(scales$conc))
  scales$conc[missing_scales_conc] <- scales_default$conc[missing_scales_conc]
#check for missing scales_time
  missing_scales_time <- setdiff(names(scales_default$time),
                                 names(scales$time))
  scales$time[missing_scales_time] <- scales_default$time[missing_scales_time]


  #fill in any non-specified optimx settings with their defaults
  optimx_settings_default = list(
    method = "bobyqa",
    itnmax = 1e6,
    hessian = FALSE,
    control = list(kkt = FALSE)
  )

  missing_optimx_settings <- setdiff(names(optimx_settings_default),
                                   names(optimx_settings))
  optimx_settings[missing_optimx_settings] <- optimx_settings_default[missing_optimx_settings]

#Create the initial pk object
  obj <- list("data_original" = data,
              "data_settings" = c(list(mapping = mapping),
                                    data_settings),
              "scales" = scales,
              "optimx_settings" = optimx_settings
              )
#nd assign it class pk
  class(obj) <- c(class(obj), "pk")
#return it
  return(obj)

}

#' Check whether an object is of class `pk`
#'
#' @param obj The object whose class is to be tested
#' @return TRUE if the object conforms to expectations for class
#'   `pk`, FALSE if it does not
#' @export
#' @author Caroline Ring
is.pk <- function(obj){
  return(inherits(obj, "pk"))
}


"+.pk" <- function(e1, e2) {
  if (missing(e2)) {
    cli::cli_abort(c(
      "Cannot use {.code +} with a single argument",
      "i" = "Did you accidentally put {.code +} on a new line?"
    ))
  }

  # Get the name of what was passed in as e2, and pass along so that it
  # can be displayed in error messages
  e2name <- deparse(substitute(e2))

  if      (is.pk(e1))  add_pk(e1, e2, e2name)
  else if (is.pkproto(e1)) {
    cli::cli_abort(c(
      "Cannot add {.cls pkproto} objects together",
      "i" = "Did you forget to add this object to a {.cls pk} object?"
    ))
  }
}

#' Add various `pkproto` objects to a `pk` object
#'
#'@param pk_obj The `pk` object
#'@param object The `pkproto` object to be added
#'@param objectname The name of the `pkproto` object to be added
#'
#'@return The `pk` object modified by the addition.
add_pk <- function(pk_obj, object, objectname) {
  if (is.null(object)) return(pk_obj)

  p <- pk_add(object, pk_obj, objectname)
  p
}

#' Add a `pk_scales` object.
#'
#' @param object The `pk_scales` object to be added.
#' @param pk_obj The `pk` object to which the `pk_scales` object will be added.
#' @param objectname The name of the `pk_scales` object.
#'
#' @return The `pk` object, modified by adding the scale.
#' @author Caroline Ring
pk_add.pk_scales <- function(object, pk_obj, objectname){
pk_obj$scales[[objectname]] <- object
}

#' Add a `pk_settings` object.
#'
#' @param object The `pk_settings` object to be added.
#' @param pk_obj The `pk` object to which the `pk_settings` object will be added.
#' @param objectname The name of the `pk_scales` object.
#'
#' @return The `pk` object, modified by adding the scale.
#' @author Caroline Ring
pk_add.pk_settings <- function(object, pk_obj, objectname){

}


pk_add.pk_stat_model <- function(object, pk_obj, objectname){

}

pk_add.pk_stat_nca <- function(object, pk_obj, objectname){

}

pk_add.pk_gof_rmse <- function(object, pk_obj, objectname){

}

pk_add.gof_r2 <- function(object, pk_obj, objectname){

}

#' Create a new `pk` object
#'
#' [pk()] initializes a new `pk` object.
#'
#'
#' [pk()] is used to construct the initial `pk` object for analysis. It is almost
#' always followed by `+` to add steps to the workflow. For example, you could
#' use `pk(my_data) + stat_model(model = '1comp')` to set up for fitting a
#' 1-compartment model.
#'
#' # The `pk` object
#'
#' A `pk` object consists of a set of concentration-dose-time data to be fitted,
#' and sets of instructions for steps in the analysis workflow:
#'
#' - settings for how to pre-process the data (harmonizing variable names, imputing missing data, calculating derived variables)
#' - scalings/transformations to be applied to the data
#' - settings for the numerical optimization algorithm to be used to fit any model
#' - optionally: which PK model(s) should be fitted to this dataset. (You do not have to fit any PK model if you don't want to; you can instead just set up the `pk` object with data, and do non-compartmental analysis on it.)
#'
#'
#' No data processing, model fitting, or any other analysis is done until you
#' explicitly request it. Until then, the `pk` object remains just a set of data
#' and instructions. This allows you to specify the instructions for each
#' analysis step without regard for the actual order of the analysis steps, and
#' to overwrite previous instructions, without having to re-do the model fitting
#' each time you add or change a set of instructions. This is particularly useful
#' if you are working in interactive mode at the R console.
#'
#'
#'
#' # Mappings
#'
#' Your input data can have any variable names you like. However, internally,
#' `invivopkfit` needs to use a set of "standard", harmonized variable names
#' (e.g., it refers to the variable containing measured tissue concentrations as
#' `Conc`; the variable containing observed time points as `Time`; and the
#' variable containing administered doses as `Dose`). In effect, `invivopkfit`
#' needs to rename the input data, and produce a new `data.frame` that uses these
#' internal harmonized variable names.
#'
#' In order to know which variable names in the input data correspond to each of
#' the internal harmonized variable names, we need to set up a mapping between
#' the internal harmonized variable names and the original variable names.
#'
#' The simplest, most flexible way to set up this mapping is by (ab)using a call
#' to [ggplot2::aes()]. In the context of [ggplot2::ggplot2-package()], you would
#' use [ggplot2::aes()] to set up mappings to `ggplot2`'s "aesthetics", internal
#' harmonized variable names which it uses for plotting: *e.g.*, `x`, `y`,
#' `color`, `size`, `shape`, and so on. In the context of
#' [invivopkfit-package()], we are setting up mappings to `invivopkfit`'s
#' internal harmonized variable names which it uses in model fitting. These
#' "`invivopkfit` aesthetic" variables are as follows:
#'
#' \itemize{
#' \item `Chemical`: A `character` variable containing the chemical identifier. All rows of `data` should have the same value for `Chemical`.
#' \item `Species`: A `character` variable containing the name of the species for which the data was measured.  All rows of `data` should have the same value for `Species`.
#' \item `Reference`: A `character` variable containing a unique identifier for the data reference (e.g., a single publication).
#' \item `Subject`: A `character` variable containing a unique identifier for the subject associated with each observation (an individual animal or group of animals).
#' \item `N_Subjects`: A `numeric` variable; an integer giving the number of individual animals represented by this observation. (Some data sets report tissue concentrations for individual animals, in which case `N_Subjects` will be 1; others report average tissue concentrations for groups of multiple animals, in which case `N_Subjects` will be greater than 1.)
#' \item `Weight`: A `numeric` variable giving the subject's body weight.
#' \item `Weight.Units`: A `character` variable giving the units of body weight.
#' \item `Route`: A `character` variable denoting the route of administration. Either `po` (oral) or `iv` (intravenous). Other routes are not currently supported.
#' \item `Dose`: A `numeric` variable giving the dose administered.
#' \item `Dose.Units`: A `character` variable giving the units of the administered doses.
#' \item `Time`: A `numeric` variable giving the time of tissue collection.
#' \item `Time.Units`: A `numeric` variable giving the units of `Time`.
#' \item `Media`: A `character` variable giving the tissue that was analyzed. Either `blood`, `plasma`, or `excreta`. Other tissues are not currently supported.
#' \item `Value`: A `numeric` variable giving the tissue concentration in units of mg/L. If `N_Subjects > 1`, `Value` is assumed to represent the mean tissue concentration for this group of subjects. If the tissue concentration was below the limit of quantification (LOQ), this value may be `NA_real_`.
#' \item `Value_SD`: A `numeric` variable giving the standard deviation of the tissue concentration in units of mg/L, if available and relevant. If `N_Subjects > 1`, `Value_SD` is assumed to represent the standard deviation of tissue concentrations for this group of subjects. If `N_Subjects == 1`, then `Value_SD` may be `NA_real_`.
#' \item `LOQ`: A `numeric` variable giving the limit of quantification applicable to this tissue concentration in units of mg/L, if available.
#' \item `Value.Units`: A `character` variable giving the units of `Value`, `Value_SD`, and `LOQ`.
#' }
#' You may additionally include mappings to other variable names of your choice,
#' which will appear in the `pk` object in `pk$data` after the analysis is done.
#'
#' As with usual calls to [ggplot2::aes()], you should provide the variable names
#' without quoting them. For example, use `ggplot2::aes(Chemical = my_chem)`. Do
#'* not* use `ggplot2::aes("Chemical" = "my_chem")`.
#'
#' Also, as with usual calls to [ggplot2::aes()], you may also specify that any
#' of the "`invivopkfit` aesthetic" variables should be mapped to a constant
#' value, rather than to a variable in `data`. For example, imagine that you
#' don't have a column in `data` that encodes the units of body weight, but you
#' know that all body weights are provided in units of kilograms. You could
#' specify `mapping = ggplot2::aes(Chemical = my_dtxsid, Species = my_species,
#' Weight = my_weight, Weight.Units = "kg")` to map `Weight.Units` to a fixed
#' value of "kg".
#'
#' Finally, as with usual calls to [ggplot2::aes()], you may specify mappings as
#' expressions that use variable names in `data`. For example, if the
#' species-name variable in `data` sometimes says "rat", sometimes "Rat",
#' sometimes "RAT", you might want to harmonize the capitalization. You can do
#' that easily by specifying `mapping = ggplot2::aes(Chemical = my_dtxsid,
#' Species = tolower(my_species)`.
#'
#' The following "aesthetics" variable names are reserved for internal use (i.e.,
#' they are automatically assigned by [preprocess_data.pk()], and should *not* be
#' included in `mapping`:
#' \itemize{
#' \item `Conc`: This is assigned as the greater of `Value` and `LOQ`, with NAs removed.
#' \item `Conc_SD`: This is set equal to `Value_SD`.
#' \item `Detect`: This is a logical variable, `TRUE` if `Conc > LOQ` and `FALSE` otherwise.
#' \item `Conc_trans`: This is `Conc` with all scalings and transformations applied as specified in `+ scale_conc()`.
#' \item `Conc_SD_trans`: This is `Conc_SD` with all scalings and transformations applied as specified in `+ scale_conc()`.
#' \item `Conc_trans.Units`: Automatically-derived from `Conc.Units` with any scalings and transformations applied. If dose normalization is requested, then `Dose.Units` is also used to automatically derive the resulting `Conc_trans.Units`. For example, if both dose-normalization and [log10()] transformation are requested, and `Conc.Units = 'mg/L'` and `Dose.Units = 'mg/kg`, then `Conc_trans.Units = log10((mg/L)/(mg/kg))`.
#' \item `Time_trans`: This is `Time` with any rescaling specified in `+ scale_time()`.
#' \item `Time_trans.Units`: The new units of time after any rescaling (e.g. `hours`, `days`, `weeks`,...)
#' }
#' If you do assign any of these reserved variable names in `mapping`, your
#' mapping will be ignored for those reserved variable names. WARNING: If you
#' have any variables with these reserved names in your original data, those
#' original variables will be dropped by [preprocess_data.pk()].
#'
#' The default value of `mapping` is the following (which refers to original
#' variable names in the built-in dataset [cvt]):
#'
#' \preformatted{
#' ggplot2::aes(
#' Chemical = chemicals_analyzed.dsstox_substance_id,
#' DTXSID = chemicals_analyzed.dsstox_substance_id,
#' CASRN = chemicals_analyzed.dsstox_casrn,
#' Species = subjects.species,
#' Reference = as.character(ifelse(is.na(documents_reference.id),
#'                                 documents_extraction.id,
#'                                 documents_reference.id)),
#' Media = series.conc_medium_normalized,
#' Route = studies.administration_route_normalized,
#' Dose = studies.dose_level_normalized,
#' Dose.Units = "mg/kg",
#' Subject = subjects.id,
#' N_Subjects =  series.n_subjects_in_series,
#' Weight = subjects.weight_kg,
#' Weight.Units = "kg",
#' Time = conc_time_values.time_hr,
#' Time.Units = "hours",
#' Value = conc_time_values.conc,
#' Value.Units = "mg/L",
#' LOQ = series.loq_normalized,
#' Value_SD  = conc_time_values.conc_sd_normalized
#' )
#' }
#'
#' # Data
#'
#' `Route` values should be either `"oral"` (oral bolus administration) or `"iv"`
#' (IV bolus administration), and `Media` values should be either `"blood"`,
#' `"plasma"`, or `"excreta"`.
#'
#' If `data` contains data for more than one `Chemical` and `Species`, then you
#' should use [facet_data()] to run a "faceted" analysis. A faceted analysis will
#' group the data according to unique combinations of the faceting variables, and
#' produce a `pk` object for each group. The result is a [tibble::tibble()]
#' grouped by the faceting variables, with a list column named `pk` containing
#' the `pk` object for each group. This [tibble::tibble()] is an object of class
#' `pk_faceted`.
#'
#' All methods for `pk` objects have a corresponding version for a
#' `pk_faceted` object, which applies the method to each `pk` object in turn and
#' either returns the same `pk_faceted` object with a modified `pk` column (for
#' methods that operate on a `pk` object and return a modified version of the
#' same `pk` object like [preprocess_data()], [data_info()], [prefit()], [fit()]), or produces a
#' [tibble::tibble()] grouped by the faceting variables, with a list column named
#' after the `pk` method containing the results of that method (for methods that
#' operate on a `pk` object but return something other than a modified `pk`
#' object, e.g. [summary.pk()], [coef.pk()], [coef_sd.pk()], [predict.pk()],
#' [residuals.pk()], [nca.pk()]).
#'
#'
#' @param data A `data.frame`. The default is an empty data frame.
#' @param mapping A mapping set up by (ab)using [ggplot2::aes()]. Call is of form
#'  `ggplot2::aes(new_variable = old_variable)` `new_variable` represents the
#'  harmonized variable name that will be used within `invivopkfit`;
#'  `old_variable` represents the variable name in `data`. If you want to
#'  provide a fixed/constant value for a `new_variable` rather than taking its
#'  value from a variable in `data`, simply supply that fixed/constant value in
#'  the `old_variable` position.
#' @param settings_preprocess_args A list of preprocessing settings.
#' @param settings_data_info_args A list of data_info settings.
#' @param settings_optimx_args A list of optimx settings.
#' @param scale_conc_args A list of concentration value scaling arguments.
#' @param scale_time_args A list of time scaling arguments
#' @param stat_model_args A list of TK model arguments.
#' @param stat_error_model_args A list of error modeling arguments
#' @param facet_data_args A list of data grouping settings.
#' @return An object of class `pk`. The initial `pk` object is a list with
#'  elements `data_orig`, `data_settings`, `scales` and `optimx_settings`.
#'  `data_orig` is the original data set to be fitted, as supplied in the
#'  argument `data`. `data_settings` is a named list containing all the other
#'  input arguments: these provide settings that will be used when the data is
#'  pre-processed before fitting.
#' @import tibble
#' @author Caroline Ring, Gilberto Padilla Mercado
#' @export

pk <- function(data = NULL,
               mapping = ggplot2::aes(
                 Chemical = analyte_dtxsid,
                 Chemical_Name = analyte_name_original,
                 DTXSID = analyte_dtxsid,
                 CASRN = analyte_casrn,
                 Species = species,
                 Reference = fk_extraction_document_id,
                 Media = conc_medium_normalized,
                 Route = administration_route_normalized,
                 Dose = invivPK_dose_level,
                 Dose.Units = "mg/kg",
                 Subject_ID = fk_subject_id,
                 Series_ID = fk_series_id,
                 Study_ID = fk_study_id,
                 ConcTime_ID = conc_time_id,
                 N_Subjects = n_subjects_normalized,
                 Weight = weight_kg,
                 Weight.Units = "kg",
                 Time = time_hr,
                 Time.Units = "hours",
                 Value = invivPK_conc,
                 Value.Units = "mg/L",
                 Value_SD = invivPK_conc_sd,
                 LOQ = invivPK_loq
               ),
               settings_preprocess_args = list(),
               settings_data_info_args = list(),
               settings_optimx_args = list(),
               scale_conc_args = list(),
               scale_time_args = list(),
               stat_model_args = list(),
               stat_error_model_args = list(),
               facet_data_args = list()

) {

  # Check to ensure the mapping contains all required harmonized column names
  mapping_default <- ggplot2::aes(
    Chemical = NA_character_,
    Species = NA_character_,
    Reference = NA_character_,
    Media = NA_character_,
    Route = NA_character_,
    Dose = NA_real_,
    Dose.Units = "mg/kg",
    Series_ID = NA_character_,
    N_Subjects = NA_real_,
    Weight = NA_real_,
    Weight.Units = "kg",
    Time = NA_real_,
    Time.Units = "hours",
    Value = NA_real_,
    Value_SD = NA_real_,
    LOQ = NA_real_,
    Value.Units = "mg/L"
  )

  missing_aes <- setdiff(names(mapping_default), names(mapping))
  if (!(length(missing_aes) == 0)) {
    mapping[missing_aes] <- mapping_default[missing_aes]
    missing_mapping <- sapply(mapping_default[missing_aes],
                              rlang::as_label)
    missing_map_names <- missing_aes
    missing_map_print <- paste(
      paste(missing_map_names, missing_mapping, sep = " = "),
      collapse = "\n"
    )
    warning(paste("'mapping' is missing the following required harmonized variables, which will be added according to the following mapping:",
                 missing_map_print,
                  sep = "\n"))
  }


  # Create the initial pk object
  obj <- list("data_original" = data,
              "mapping" = mapping,
              "status" = status_init
  )

  # nd assign it class pk
  class(obj) <- append("pk", class(obj))

  # Add default data preprocessing settings
  obj <- obj + do.call(settings_preprocess, settings_preprocess_args)

  # Add default data info settings
  obj <- obj + do.call(settings_data_info,
                       settings_data_info_args)

  # Add default optimx settings
  obj <- obj + do.call(settings_optimx,
                       settings_optimx_args)

  # Add default scalings for conc and time
  obj <- obj + do.call(scale_conc, scale_conc_args)
  obj <- obj + do.call(scale_time, scale_time_args)

  # #Add default models: flat, 1comp, 2comp
  obj <- obj + do.call(stat_model, stat_model_args)

  # Add default error model
  obj <- obj + do.call(stat_error_model, stat_error_model_args)

  # add default faceting
  obj <- obj + do.call(facet_data, facet_data_args)

  # return the initialized pk object
  return(obj)
}

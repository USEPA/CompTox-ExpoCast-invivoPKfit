#' Main fitting function
#'
#' Fits parameters of a specified model to concentration-time data.
#'
#' @param data.set A \code{data.frame} of concentration-time data. Preferably
#'   \code{pkdataset_nheerlcleaned} or \code{pkdataset_nheerlorig}.
#' @param model Which general model should be fit for each chemical. Presently,
#'   only `"flat"`, `"1compartment"`, and `"2compartment"` are implemented.
#' @param modelfun Character, either `"analytic"` or `"full"` -- whether to fit
#'   using the analytic solution to the model, or numerical solution of the full
#'   ODE model. Presently, `"analytic"` is the default, and is recommended
#'   because the analytic solutions are exact and much faster than numerically
#'   integrating the full ODE model.
#' @param preprocess Logical: `TRUE` to preprocess data using
#'   [preprocess_data()]; `FALSE` to use already-pre-processed data (i.e. the
#'   output of [preprocess_data()]. Default `TRUE`.
#' @param names_list Used only if `preprocess == TRUE`. As for
#'   [preprocess_data()] and [rename_columns()]: A named list where the names
#'   are the new variable names, and the values are the old variable names. If a
#'   value is NULL, then there is no old variable corresponding to that new
#'   variable. A new variable will be added, with default value as given in
#'   `defaults_list`; if no default is given in `defaults_list`, the new
#'   variable will be filled with `NA_character`.
#' @param defaults_list Used only if `preprocess == TRUE`. As for
#'   [preprocess_data()] and [rename_columns()]: A named list where the names
#'   are the new variable names, and the values are the default values to fill
#'   for each new variable, if the corresponding old variable name is NULL in
#'   `names_list`.
#' @param ratio_conc_to_dose Used only if `preprocess == TRUE`. As for
#'   [preprocess_data()]. Ratio between the mass units used to report the
#'   concentration data and the mass units used to report the dose. Default 1.
#'   For example, concentration reported in ug/L and dose reported in mg/kg/day
#'   would require `ratio_conc_to_dose = 0.001`, because 1 ug/1 mg = 1e-6 g /
#'   1e-3 g = 0.001.
#' @param calc_loq_factor Used only if `preprocess == TRUE`. As for
#'   [preprocess_data()] and [estimate_loq()].
#' @param routes_keep Used only if `preprocess == TRUE`. As for
#'   [preprocess_data()]. Character vector of administration routes to keep.
#'   Default `c("po", "iv")` to keep only oral and intravenous data (because at
#'   present, model fitting is only implemented for oral and IV routes). Any
#'   data where route is not on this list will be discarded. (For example, with
#'   the default, inhalation and dermal data will be discarded.)
#' @param media_keep Used only if `preprocess == TRUE`. As for
#'   [preprocess_data()]. Character vector of sample media to keep. Default is
#'   `c("blood", "plasma")` to keep only concentration data measured in blood
#'   and plasma (because at present, model fitting is only implemented for blood
#'   and plasma). Any data where medium is not on this list will be discarded.
#'   (For example, with the default, urine data will be discarded.)
#' @param impute_loq Used only if `preprocess == TRUE`. As for
#'   [preprocess_data()]. Logical: TRUE to impute values for missing LOQs; FALSE
#'   to leave them alone.
#' @param impute_sd Used only if `preprocess == TRUE`. As for
#'   [preprocess_data()]. Logical: TRUE to impute values for missing sample SDs
#'   for multi-subject observations; FALSE to leave them alone
#' @param study_def Used only if `preprocess == TRUE`. As for
#'   [preprocess_data()]. A character vector specifying the variables (new names)
#'   whose unique combinations define individual "studies" (where each "study"
#'   will have its own error SD in the fitting process). Default is `c("DTXSID",
#'   "Species", "Reference", "Route", "Media")`.
#'@param pool_sigma Logical: Whether to pool all data (estimate only one error
#'  standard deviation) or not (estimate separate error standard deviations for
#'  each study). Default FALSE to estimate separate error SDs for each study.
#'  (If `fitdata` only includes one study, `pool_sigma` will have no effect,
#'  because only one error SD would be estimated in the first place.)
#'@param fit_conc_dose Logical: Whether to fit dose-normalized concentrations
#'  (TRUE) or non-dose-normalized concentrations (FALSE). Default TRUE.
#' @param rescale_time As for `analyze_subset()`. Logical: TRUE to rescale time
#'   from hours if latest detection time corresponds to days, weeks, months, or
#'   years. FALSE to leave time in units of hours. Default TRUE.
#' @param get_starts_args Named list or NULL: any additional arguments to
#'   [get_starts()] (other than `model` and `fitdata`, which are always passed).
#'   Default NULL to accept the default arguments for [get_starts()].
#' @param get_lower_args Named list or NULL: any additional arguments to
#'   [get_lower_bounds()] (other than `model` and `fitdata`, which are always
#'   passed). Default NULL to accept the default arguments for
#'   [get_lower_bounds()].
#' @param get_upper_args Named list or NULL: any additional arguments to
#'   [get_upper_bounds()] (other than `model` and `fitdata`, which are always
#'   passed). Default NULL to accept the default arguments for
#'   [get_upper_bounds()].
#' @param optimx_args As for [analyze_subset()]. A named list of additional
#'   arguments to [optimx::optimx()] (used for parameter estimation), other than
#'   `par`, `fn`, `lower`, and `upper`. Default is:
#'
#'    ```
#'     list(
#'           "method" = "bobyqa",
#'           "itnmax" = 1e6,
#'           "control" = list("kkt" = FALSE)
#'          )
#'    ```
#'   Briefly:  `"method"` allows you to select an optimization algorithm.
#'   `"itnmax"` controls the maximum number of iterations allowed in an attempt
#'   to optimize. `"control"` is itself a named list of control parameters
#'   relevant to the selected method. See documentation for [optimx::optimx()]
#'   for more details and more options. Note lower and upper bounds (box
#'   constraints) will be supplied; if you want them to be respected, please
#'   choose a method that allows box constraints (e.g. "bobyqa" or "L-BFGS-B").
#' @param suppress.messages Logical: Whether to suppress verbose messages.
#'   Default FALSE, to be verbose.
#'
#' @return A `data.table` of fitted parameter values for each dataset.
#'
#' @author Caroline Ring, John Wambaugh
#' @export fit_all
#' @import data.table

fit_all <- function(data.set,
                    model,
                    modelfun = "analytic",

                    preprocess = TRUE,

                    names_list =list(
                      "Compound_Dosed" = "studies.test_substance_name_original",
                      "DTXSID_Dosed" = "chemicals_dosed.dsstox_substance_id",
                      "CAS_Dosed" = "chemicals_dosed.dsstox_casrn",
                      "Compound_Analyzed" = "series.analyte_name_original",
                      "DTXSID_Analyzed" = "chemicals_analyzed.dsstox_substance_id",
                      "CAS_Analyzed" = "chemicals_analyzed.dsstox_casrn",
                      "Reference" = "documents_reference.id",
                      "Extraction" = "documents_extraction.id",
                      "Species" = "subjects.species",
                      "Weight" ="subjects.weight_kg",
                      "Weight.Units" = NULL,
                      "Dose" = "studies.dose_level_normalized",
                      "Dose.Units" = NULL,
                      "Time" = "conc_time_values.time_hr",
                      "Time.Units" = NULL,
                      "Media" = "series.conc_medium_normalized",
                      "Value" = "conc_time_values.conc",
                      "Value.Units" = NULL,
                      "Route" = "studies.administration_route_normalized",
                      "LOQ" = "series.loq",
                      "Subject" = "subjects.id",
                      "N_Subjects" = "series.n_subjects_in_series",
                      "Value_SD" = "conc_time_values.conc_sd"),
                    defaults_list =   list(
                      "Weight.Units" = "kg",
                      "Dose.Units" = "mg/kg",
                      "Time.Units" = "hours",
                      "Value.Units" = "mg/L"),

                    ratio_conc_to_dose = 1,
                    calc_loq_factor = 0.45,
                    routes_keep = c("po", "iv"),
                    media_keep = c("blood", "plasma"),
                    impute_loq = TRUE,
                    impute_sd = TRUE,
                    study_def = c("DTXSID", "Species", "Reference", "Route", "Media"),

                    fit_conc_dose = TRUE,
                    fit_log_conc = FALSE,
                    rescale_time = TRUE,
                    get_starts_args = list(start_from_httk = "all",
                                           start_from_data = "all"),
                    get_lower_args = list(Vdist_from_species = FALSE),
                    get_upper_args = list(Fgutabs_Vdist_from_species = FALSE,
                                          sigma_from_data = TRUE),
                    optimx_args = list(
                      "method" = "bobyqa",
                      "itnmax" = 1e6,
                      "control" = list("kkt" = FALSE)
                    ),

                    suppress.messages = FALSE) {

  #############################
  #############################
  #############################
  ## PREPROCESS DATA #########
  #############################
  #############################
  #############################

  if(preprocess %in% TRUE){
 data.set <- preprocess_data(data.set = data.set,
                             names_list = names_list,
                             defaults_list = defaults_list,
                             ratio_conc_to_dose = ratio_conc_to_dose,
                             calc_loq_factor = calc_loq_factor,
                             routes_keep = routes_keep,
                             media_keep = media_keep,
                             impute_loq = impute_loq,
                             impute_sd = impute_sd,
                             study_def = study_def,
                             suppress.messages = suppress.messages)
  }

 data.set <- as.data.table(data.set)

  #Non-comapartmental fits:
  if (model=="noncompartment") {

    PK.fit.table <- do_noncomp_fit(data.set)

  } else {

    #Analyze by chemical & species first.
if(!suppress.messages){
  message("Analyzing data by chemical and species, with each study having its own error SD...")
}
    PK.fit.joint <- data.set[,
                             analyze_subset(fitdata = .SD,
                                            #fitdata = .SD passes a data.table
                                            #consisting of the columns
                                            #referenced in .SDcols, split by the
                                            #values in "by"
                                             modelfun = modelfun,
                                             model = model,
                                            pool_sigma = FALSE,
                                            fit_conc_dose = fit_conc_dose,
                                            fit_log_conc = fit_log_conc,
                                            rescale_time = rescale_time,
                                            get_starts_args = get_starts_args,
                                            get_lower_args = get_lower_args,
                                            get_upper_args = get_upper_args,
                                            optimx_args = optimx_args,
                                             suppress.messages=suppress.messages),
                             .SDcols = names(data.set),
                             #.SDcols argument: By default, .SD includes all
                             #columns *except* those named in "by=". We want to
                             #*include* those named in "by=", so we explicitly
                             #state that .SDcols must include all columns in
                             #data.set. Otherwise, "fitdata" would be missing
                             #the "DTXSID" and "Species" columns.
                             by = c("DTXSID", #split data into unique combinations of these columns
                                    "Species")
                             ]

    ### If chemical/species combination has multiple studys:
    ### then analyze each study individually,
    ###and also do a pooled analysis without individual study sigmas
    data.set[, MultipleStudies := length(unique(Study)) > 1,
             by = c("DTXSID", "Species")]
    multi.study.cas <- unique(subset(data.set,
                                   MultipleStudies == TRUE)$DTXSID)
    if (length(multi.study.cas) > 0) {
      if(!suppress.messages){
        message("Analyzing data by chemical, species, and study...")
      }
      data.set.multi.study <- subset(data.set, MultipleStudies == TRUE)

      PK.fit.separate <- data.set.multi.study[,
                                 analyze_subset(fitdata = .SD,
                                                 modelfun = modelfun,
                                                 model = model,
                                                pool_sigma = FALSE,
                                                fit_conc_dose = fit_conc_dose,
                                                fit_log_conc = fit_log_conc,
                                                rescale_time = rescale_time,
                                                get_starts_args = get_starts_args,
                                                get_lower_args = get_lower_args,
                                                get_upper_args = get_upper_args,
                                                optimx_args = optimx_args,
                                                suppress.messages = suppress.messages),
                                 .SDcols = names(data.set.multi.study),
                                                 by = c(
                                                   "DTXSID",
                                                   "Species",
                                                   "Study")]

      PK.fit.separate[, Study := NULL]


      if(!suppress.messages){
        message("Analyzing data by chemical and species, pooling all studies...")
      }
      PK.fit.pooled <- data.set.multi.study[,
                                          analyze_subset(fitdata = .SD,
                                                         modelfun = modelfun,
                                                         model = model,
                                                         pool_sigma = TRUE,
                                                         rescale_time = rescale_time,
                                                         fit_conc_dose = fit_conc_dose,
                                                         fit_log_conc = fit_log_conc,
                                                         get_starts_args = get_starts_args,
                                                         get_lower_args = get_lower_args,
                                                         get_upper_args = get_upper_args,
                                                         optimx_args = optimx_args,
                                                         suppress.messages = suppress.messages),
                                          .SDcols = names(data.set.multi.study),
                                          by = c(
                                            "DTXSID",
                                            "Species")]

      PK.fit.bind <- rbindlist(list("Joint" = PK.fit.joint,
                           "Separate" = PK.fit.separate,
                            "Pooled" = PK.fit.pooled
      ),
                           use.names = TRUE,
                           fill = TRUE,
                           idcol = "Analysis_Type")
    } else {
      PK.fit.bind <- rbindlist(list("Joint" = PK.fit.joint
      ),
      use.names = TRUE,
      fill = TRUE,
      idcol = "Analysis_Type")
    }

    #record which model was fit to data and whether it was full or analytical
    PK.fit.bind[, model := model]
    PK.fit.bind[, model.type := modelfun]

    #############################
    #############################
    #############################
}

  return(PK.fit.bind)
}

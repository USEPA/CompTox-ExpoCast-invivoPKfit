#' Main fitting function
#'
#' Fits parameters of a specified model to concentration-time data.
#'
#' @param data.set A \code{data.frame} of concentration-time data. Preferably
#'   \code{pkdataset_nheerlcleaned} or \code{pkdataset_nheerlorig}.
#' @param model Which general model should be fit for each chemical. Presently,
#'   only "1compartment" and "2compartment" are implemented.
#' @param modelfun Either "analytic" or "full" -- whether to fit using the
#'   analytic solution to the model, or the full ODE model. Presently,
#'   "analytic" is recommended (because the analytic solution is exact and much
#'   faster).
#' @param ratio_conc_to_dose Ratio between the mass units used to report the
#'   concentration data and the mass units used to report the dose. Default 1.
#'   For example, concentration reported in ug/L and dose reported in mg/kg/day
#'   would require \code{ratio_conc_to_dose = 0.001}, because 1 ug/1 mg = 1e-6 g
#'   / 1e-3 g = 0.001.
#' @param names_list As for [rename_columns()]: A named list where the names are
#'   the new variable names, and the values are the old variable names. If a
#'   value is NULL, then there is no old variable corresponding to that new
#'   variable. A new variable will be added, with default value as given in
#'   `defaults_list`; if no default is given in `defaults_list`, the new
#'   variable will be filled with `NA_character`.
#' @param defaults_list As for [rename_columns()]: A named list where the names
#'   are the new variable names, and the values are the default values to fill
#'   for each new variable, if the corresponding old variable name is NULL in
#'   `names_list`.
#' @param get_starts_args Named list or NULL: any additional arguments to [get_starts()] (other than
#'   `model` and `fitdata`, which are always passed). Default NULL to accept the
#'   default arguments for [get_starts()].
#' @param get_lower_args Named list or NULL: any additional arguments to [get_lower_bounds()] (other
#'   than `model` and `fitdata`, which are always passed). Default NULL to
#'   accept the default arguments for [get_lower_bounds()].
#' @param get_upper_args Named list or NULL: any additional arguments to [get_upper_bounds()] (other
#'   than `model` and `fitdata`, which are always passed). Default NULL to
#'   accept the default arguments for [get_upper_bounds()].
#' @param optimx_args A named list of additional arguments to [optimx::optimx()]
#'   (used for parameter estimation), other than `par`, `fn`, `lower`, and
#'   `upper`. Default is:
#'
#'    ```
#'     list(
#'           "method" = "bobyqa",
#'           "itnmax" = 1e6,
#'           "control" = list("maximize" = TRUE,
#'                            "kkt" = FALSE)
#'          )
#'    ```
#'  See documentation for [optimx::optimx()] for possible arguments and details.
#'  Note lower and upper bounds (box constraints) will be supplied; if you want
#'  them to be respected, please choose a method that allows box constraints
#'  (e.g. "bobyqa" or "L-BFGS-B").
#' @param suppress.messages Logical: Whether to suppress verbose messages.
#'   Default FALSE, to be verbose.
#'
#' @return A `data.table` of fitted parameter values for each chemical.
#'
#' @author Caroline Ring, John Wambaugh
#' @export fit_all
#' @importFrom PK nca.batch
#' @importFrom magrittr "%>%"

fit_all <- function(data.set,
                    model,
                    modelfun = NA,

                    names_list =list(
                      "Compound" = "chemicals_dosed.id",
                      "DTXSID" = "chemicals_dosed.dsstox_substance_id",
                      "CAS" = "chemicals_dosed.dsstox_casrn",
                      "Reference" = "documents_reference.id",
                      "Extraction" = "documents_extraction.id",
                      "Weight" ="subjects.weight_kg",
                      "Weight.Units" = NULL,
                      "Dose" = "studies.dose_level_normalized",
                      "Dose.Units" = NULL,
                      "Time" = "conc_time_values.time_hr",
                      "Time.Units" = NULL,
                      "Media" = "series.conc_medium_normalized",
                      "Value" = "series.conc",
                      "Value.Units" = NULL,
                      "Route" = "studies.administration_route_normalized",
                      "Extraction" = "documents_extraction.id",
                      "LOQ" = "series.loq",
                      "Subject" = "subjects.id"),
                    defaults_list =   list(
                      "Weight.Units" = "kg",
                      "Dose.Units" = "mg/kg",
                      "Time.Units" = "h",
                      "Value.Units" = "mg/L"),

                    ratio_conc_to_dose = 1,
                    calc_loq_factor = 0.45,

                    get_starts_args = NULL,
                    get_lower_args = NULL,
                    get_upper_args = NULL,
                    optimx_args = list(
                      "method" = "bobyqa",
                      "itnmax" = 1e6,
                      "control" = list("maximize" = TRUE,
                                       "kkt" = FALSE)
                    ),

                    suppress.messages = FALSE) {

  #############################
  #############################
  #############################
  ## PREPROCESS DATA #########
  #############################
  #############################
  #############################

 data.set <- preprocess_data(data.set = data.set,
                             names_list = names_list,
                             defaults_list = defaults_list,
                             ratio_conc_to_dose = ratio_conc_to_dose,
                             calc_loq_factor = calc_loq_factor,
                             suppress.messages = suppress.messages)

  #Non-comapartmental fits:
  if (model=="noncompartment") {

    PK.fit.table <- do_noncomp_fit(data.set)

  } else {

    #Analyze by chemical & species first.
if(!suppress.messages){
  message("Analyzing data by chemical and species, with each reference having its own error SD...")
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

    ### If chemical/species combination has multiple references:
    ### then analyze each reference individually,
    ###and also do a pooled analysis without individual reference sigmas
    data.set[, MultipleReferences := length(unique(Reference)) > 1,
             by = c("DTXSID", "Species")]
    multi.ref.cas <- unique(subset(data.set,
                                   MultipleReferences == TRUE)$DTXSID)
    if (length(multi.ref.cas) > 0) {
      if(!suppress.messages){
        message("Analyzing data by chemical, species, and reference...")
      }
      data.set.multi.ref <- subset(data.set, MultipleReferences == TRUE)
      data.set.multi.ref[, Reference_orig := Reference]
      PK.fit.separate <- data.set.multi.ref[,
                                 analyze_subset(fitdata = .SD,
                                                 modelfun = modelfun,
                                                 model = model,
                                                pool_sigma = FALSE,
                                                get_starts_args = get_starts_args,
                                                get_lower_args = get_lower_args,
                                                get_upper_args = get_upper_args,
                                                optimx_args = optimx_args,
                                                suppress.messages = suppress.messages),
                                 .SDcols = names(data.set.multi.ref),
                                                 by = c(
                                                   "DTXSID",
                                                   "Species",
                                                   "Reference_orig")]


      if(!suppress.messages){
        message("Analyzing data by chemical and species, pooling all references...")
      }
      PK.fit.pooled <- data.set.multi.ref[,
                                          analyze_subset(fitdata = .SD,
                                                         modelfun = modelfun,
                                                         model = model,
                                                         pool_sigma = TRUE,
                                                         get_starts_args = get_starts_args,
                                                         get_lower_args = get_lower_args,
                                                         get_upper_args = get_upper_args,
                                                         optimx_args = optimx_args,
                                                         suppress.messages = suppress.messages),
                                          .SDcols = names(data.set.multi.ref),
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
    } else PK.fit.bind <- PK.fit.joint

    #record which model was fit to data and whether it was full or analytical
    PK.fit.bind[, model := model]
    PK.fit.bind[, model.type := modelfun]

    #############################
    #############################
    #############################
}

  return(PK.fit.bind)
}

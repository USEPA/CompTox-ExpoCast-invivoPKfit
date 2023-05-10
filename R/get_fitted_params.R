#' Get a set of fitted model parameters
#'
#' Get a set of fitted model parameters for a specific dataset, analysis, and
#' model
#'
#' [fit_all()] produces a table of fitted parameters for all chemicals and
#' species, for a given analysis and model. [merge_fits()] post-processes and
#' merges these tables for multiple models and fitting conditions. This
#' facilitates tabular display of fitted parameters. However, it can also be
#' useful to pull a particular set of fitted parameters from the table, in a
#' form suitable for input to one of the model functions, e.g. [cp_flat()],
#' [cp_1comp()], or [cp_2comp()]. This utility function does the necessary data
#' wrangling.
#'
#' The user provides a full, merged table of PK parameters as produced by
#' [merge_fits()]. The user also specifies which fit they want: a unique
#' combination of model, fitting conditions (log-transformation,
#' dose-normalization, rescaling time, optimization method used), substance
#' (DTXSID), species, analysis type (Joint, Separate, or Pooled), and the list
#' of studies analyzed (only necessary for Separate analyses, since these were
#' done one study at a time).
#'
#' @param pk_fit A merged table of fitted PK parameters, as produced by [merge_fits()].
#'   Must contain variables `DTXSID`, `Species`, `Analysis_Type`, and
#'   `Studies.Analyzed`. Must also contain variables corresponding to the
#'   parameters of the model specified in argument `model_in`, as required for
#'   the route and media of data in `newdata` (as given by
#'   [get_model_paramnames()]). These variables must be named as
#'   `[param].[model]`. For example, for the 1-compartment model, if `newdata`
#'   contains only IV-dosing data measured in plasma, then `pk_fit` must contain
#'   variables named `kelim.1compartment` and `Vdist.1compartment`.
#' @param model_in The model to evaluate: one of 'flat', '1compartment', or
#'   '2compartment'.
#' @param fit_log_conc_in Whether the fit was done to log-transformed
#'   concentrations: TRUE or FALSE.
#' @param fit_conc_dose_in Whether the fit was done to dose-normalized
#'   concentrations: TRUE or FALSE.
#' @param rescale_time_in Whether the fit was done after rescaling time from
#'   hours to days, weeks, months, or years, for long studies: TRUE or FALSE.
#' @param method The algorithm used to do the fitting: "bobyqa" or "L-BFGS-B."
#' @param DTXSID_in The DSSTox Substance ID for which to evaluate
#' @param Species_in The species for which to evaluate
#' @param Analysis_Type_in The analysis type that produced the fit being
#'   evaluated: one of 'Joint', 'Separate', or 'Pooled'
#' @param Studies.Analyzed_in The comma-separated string of studies included in
#'   the fit being analyzed.
#' @return A named list of model parameters for the specified fit, suitable for
#'   input into the `params` argument of one of the model functions.
#' @author Caroline Ring
#' @export
#'
get_fitted_params <- function(pk_fit_row,
                              model_in){
  if(!(grepl(x = model_in,
             pattern = "None"))){

    #get parameter names for this model
    param_names <- get_model_paramnames(model = model_in)
    #get names of appropriate columns of pk_fit_row -- named [param].[model]
    param_names_dot <- paste(param_names, model_in, sep = ".")
    #extract appropriate columns of pk_fit_row
    params <- pk_fit_row[, .SD, .SDcols = intersect(param_names_dot,
                                                names(pk_fit_row))]
    #rename to strip the ".[model]" part
    setnames(params,
             param_names_dot,
             param_names,
             skip_absent = TRUE)
    #set Rblood2plasma to 1, if it is there and NA
    if("Rblood2plasma" %in% names(params)){
      if(is.na(params$Rblood2plasma)){
        params$Rblood2plasma <- 1
      }
    }
    #convert to list
    params <- as.list(params)
    #remove any NA params
    params <- params[!(sapply(params,
                              is.na))]

    return(params)
  }else{
    return(NULL)
  }
}

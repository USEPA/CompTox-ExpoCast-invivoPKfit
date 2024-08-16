#' Evaluate TK statistics
#'
#' Evaluate TK statistics from a fitted model by comparing to NCA results
#'
#' @param obj A [pk()] model object. Must be fitted, or the function will exit
#'   with an error.
#' @param newdata Optional: A `data.frame` containing new data for which to
#'   compute the TK stats. Must contain at least variables `Chemical`,
#'   `Species`, `Route`, `Media`, `Dose`, `Dose.Units`, `Conc.Units`, either
#'   `Time_trans.Units` or `Time.Units`, and any other variables named in
#'   `tk_grouping`. Default `NULL`, to use the data in `obj$data`.
#' @param tk_group A list of variables provided using a `dplyr::vars()` call.
#'   The data (either `newdata` or `obj$data`) will be grouped according to the
#'   unique combinations of these variables. For each unique combination of
#'   these variables in the data, a set of TK statistics will be computed. The
#'   default is `obj$settings_data_info$nca_group`, to derive TK statistics for
#'   the same groups of data as non-compartmental analysis statistics. With the
#'   default, you can directly compare e.g. a model-predicted AUC_inf to the
#'   corresponding NCA-estimated AUC_inf. However, you may specify a different
#'   data grouping if you wish. Each group should have a unique combination of
#'   `Chemical`, `Species`, `Route`, `Media`, and `Dose`, because the TK stats
#'   depend on these values, and it is required to have one unique set of TK
#'   stats per group.
#' @param model Character: One or more of the models fitted. Default `NULL` to
#'   return TK stats for all models.
#' @param method Character: One or more of the [optimx::optimx()] methods used.
#'   Default `NULL` to return TK stats for all methods.
#' @param exclude Logical: `TRUE` to get the TK groupings after removing any
#'   observations in the data marked for exclusion (if there is a variable
#'   `exclude` in the data, an observation is marked for exclusion when `TRUE`).
#'   `FALSE` to include all observations when getting the TK
#'   groupings, regardless of exclusion status. Default `TRUE`.
#' @param dose_norm Logical: `TRUE` (default) to dose-normalize before
#'   calculating both the NCA statistics and the fitted TK statistics (i.e. all
#'   dose-dependent statistics will be for a unit dose of 1 mg/kg, including
#'   Cmax, AUC, Css). `FALSE` to calculate NCA and fitted TK stats separately
#'   for each dose group (you must specify `Dose` as one of the variables in
#'   `tk_group` for this to work). If `dose_norm` is `TRUE` and you also specify
#'   `Dose` as one of the `tk_group` variables, then the dose part of the
#'   grouping will be ignored in the output. (If `TRUE`, under
#'   the hood, this function will temporarily overwrite the `Dose` column in its
#'   local copy of `newdata` with 1's. This doesn't affect the data outside of
#'   this function. But it means that any values in the `Dose` variable of
#'   `newdata` will be ignored if `TRUE`.)
#' @param finite_only Logical: `TRUE` (default) returns only rows (observations)
#'   for which AUC is finite in both `nca` and `tkstats`. This also means it will
#'   by default never return instances where winning model == `model_flat`.
#' @param ... Additional arguments. Currently not in use.
#' @return A `data.frame` with one  row for each "winning" model in
#'   `model` from [get_winning_model()]. The `data.frame` will have the variables
#'   returned by the `tkstats_fun` for its corresponding model. (For the
#'   built-in models `model_flat`, `model_1comp`, and `model_2comp`, these
#'   variables are `param_name` and `param_value`.) Additionally, there will be
#'   a variable `method` denoting the [optimx::optimx()] method used to optimize
#'   the set of model parameters used to derive each set of TK statistics.
#' @import tidyselect
#' @export
#' @author Caroline Ring, Gilberto Padilla Mercado, John Wambaugh
#' @family methods for fitted pk objects

eval_tkstats.pk <- function(obj,
                            newdata = NULL,
                            model = NULL,
                            method = NULL,
                            tk_group = NULL,
                            exclude = TRUE,
                            dose_norm = TRUE,
                            finite_only = TRUE,
                            ...){

  #ensure that the model has been fitted
  check <- check_required_status(obj = obj,
                                 required_status = status_fit)
  if(!(check %in% TRUE)){
    stop(attr(check, "msg"))
  }

  if(is.null(model)) model <- names(obj$stat_model)
  if(is.null(method)) method <- obj$settings_optimx$method
  if(is.null(newdata)) newdata <- obj$data
  if(is.null(tk_group)) tk_group <- obj$settings_data_info$nca_group

  method_ok <- check_method(obj = obj, method = method)
  model_ok <- check_model(obj = obj, model = model)

  # Grouping variables
  grp_vars <- sapply(tk_group, rlang::as_label)
  data_grp_vars <- sapply(obj$data_group, rlang::as_label)
  error_grp_vars <- sapply(obj$stat_error_model$error_group, rlang::as_label)

  newdata_ok <- check_newdata(newdata = newdata,
                              olddata = obj$data,
                              req_vars = union(
                                c("Chemical",
                                  "Species",
                                  "Time",
                                  "Time.Units",
                                  "Dose",
                                  "Conc",
                                  "Dose.Units",
                                  "Conc.Units",
                                  "Route",
                                  "Media"),
                                grp_vars),
                              exclude = exclude)

  #if exclude = TRUE, remove excluded observations
  if (exclude %in% TRUE) {
    newdata <- subset(newdata, exclude %in% FALSE)
  }

  if (dose_norm %in% TRUE) {
    newdata$Conc <- newdata$Conc/newdata$Dose
    newdata$Dose <- newdata$Dose/newdata$Dose
  } # Need this transformation for get_tkstats

  #Get the winning model for filtering
  winmodel_df <- get_winning_model(obj = obj,
                                   method = method)

  #calc NCA for newdata
  nca_df <- get_nca(obj = obj)
  nca_df <- nca_df %>% dplyr::select(-param_sd_z, -param_units) %>%
    tidyr::pivot_wider(names_from = param_name,
                       values_from = param_value) %>%
    dplyr::right_join(winmodel_df) %>%
    dplyr::relocate(model, method, .after = Media)



  #get tkstats
  tkstats_df <- get_tkstats(obj = obj,
                            newdata = newdata,
                            model = model,
                            method = method,
                            tk_group = tk_group,
                            dose_norm = dose_norm,
                            exclude = exclude) %>%
    ungroup()

  tkstats_df <- dplyr::left_join(winmodel_df %>% dplyr::ungroup(),
                                 tkstats_df,
                                 by = c(data_grp_vars, "method", "model"))

  # Assess group variables
  # Are all error_group variables in data_group, error group might be subset
  if(!all(error_grp_vars %in% data_grp_vars)) {
    remaining_vars <- error_grp_vars[!(error_grp_vars %in% data_grp_vars)]

    if (dose_norm) {
      tkstats_df_refs <- tkstats_df %>%
        dplyr::group_by(!!!obj$data_group, Route, Media) %>%
        dplyr::summarize(dplyr::across(all_of(remaining_vars),
                                       \(x) {paste0(x, collapse = ",")}))
      nca_df <- nca_df %>%
        dplyr::ungroup() %>%
        dplyr::select(!c(Dose, design))
    } else {
      tkstats_df_refs <- tkstats_df %>%
        dplyr::group_by(!!!obj$data_group, Route, Media, Dose) %>%
        dplyr::summarize(dplyr::across(all_of(remaining_vars),
                                       \(x) {paste0(x, collapse = ",")}))
    }


    tkstats_df <- dplyr::left_join(tkstats_df_refs,
                                   tkstats_df %>%
                                     dplyr::select(!all_of(remaining_vars))) %>%
      dplyr::ungroup() %>%
      dplyr::distinct()
  }


  # prepare for merge
  nca_df_red <- nca_df %>%
    dplyr::rename_with(~ paste0(.x, ".nca", recycle0 = TRUE),
                       !any_of(c(grp_vars, "model", "method")))

  tkstats_df_red <- tkstats_df %>%
    dplyr::rename_with(~ paste0(.x, ".tkstats", recycle0 = TRUE),
                       !any_of(c(grp_vars, "Dose.Units", "Conc.Units",
                                 "Time.Units", "Time_trans.Units",
                                 "Rblood2plasma", "model", "method")))

  if (dose_norm) {
    message("TK statistics calculated based on dose-normalized value of 1mg/kg")
  }

  # Merge the tkstats and the nca data.frames
    tk_eval <- suppressMessages(dplyr::left_join(tkstats_df_red,
                                nca_df_red))

  # Filter out infinite values for AUC_infinity
  if (finite_only) {
    tk_eval <- tk_eval %>%
      filter(is.finite(AUC_infinity.tkstats),
             is.finite(AUC_infinity.nca))
  }


  return(tk_eval)


}

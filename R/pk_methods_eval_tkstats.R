#' Evaluate TK statistics
#'
#' Evaluate TK statistics from a fitted model by comparing to NCA results
#'
#' @inheritParams get_tkstats.pk
#' @param finite_only Logical (Default: TRUE). If FALSE, will include non-finite values
#'   for `AUC_infinity` from both compartmental and noncompartmental analysis.
#' @return A `data.frame` with one  row for each "winning" model in
#'   `model` from [get_winning_model()]. The `data.frame` will have the variables
#'   returned by the `tkstats_fun` for its corresponding model. (For the
#'   built-in models `model_flat`, `model_1comp`, and `model_2comp`, these
#'   variables are `param_name` and `param_value`.) Additionally, there will be
#'   a variable `method` denoting the [optimx::optimx()] method used to optimize
#'   the set of model parameters used to derive each set of TK statistics.
#' @import dplyr
#' @export
#' @author Caroline Ring, Gilberto Padilla Mercado, John Wambaugh
#' @family methods for fitted pk objects

eval_tkstats.pk <- function(obj,
                            newdata = NULL,
                            model = "winning",
                            method = NULL,
                            tk_group = NULL,
                            exclude = TRUE,
                            dose_norm = FALSE,
                            finite_only = FALSE,
                            suppress.messages = NULL,
                            ...) {

  if (is.null(suppress.messages)) {
    suppress.messages <- obj$pk_settings$preprocess$suppress.messages
  }

  # ensure that the model has been fitted
  check <- check_required_status(obj = obj, required_status = status_fit)

  if (!(check %in% TRUE)) {
    stop(attr(check, "msg"))
  }

  if (is.null(model)) model <- names(obj$stat_model)
  if (is.null(method)) method <- obj$pk_settings$optimx$method
  if (is.null(newdata)) newdata <- get_data.pk(obj)
  if (is.null(tk_group)) tk_group <- get_nca_group(obj)

  method_ok <- check_method(obj = obj, method = method)
  if (all(model %in% "winning")) {
    win <- TRUE
    model <- names(obj$stat_model)
  } else {
    win <- FALSE
  }
  model_ok <- check_model(obj = obj, model = model)

  # Grouping variables
  grp_vars <- sapply(tk_group, rlang::as_label)
  data_grp_vars <- get_data_group.pk(obj, as_character = TRUE)
  error_grp_vars <- get_error_group.pk(obj, as_character = TRUE)

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

  # if exclude = TRUE, remove excluded observations
  if (exclude %in% TRUE) {
    newdata <- subset(newdata, exclude %in% FALSE)
  }

  if (dose_norm %in% TRUE) {
    newdata$Conc <- newdata$Conc / newdata$Dose
    newdata$Dose <- newdata$Dose / newdata$Dose
  } # Need this transformation for get_tkstats

  if (win %in% TRUE) {
  # Get the winning model for filtering
    winmodel_df <- get_winning_model.pk(obj = obj,
                                        method = method) |>
      dplyr::select(-c(near_flat, preds_below_loq))
  }


  # calc NCA for newdata
  nca_df <- nca(obj = obj,
                newdata = newdata,
                nca_group = tk_group,
                dose_norm = dose_norm,
                exclude = exclude,
                suppress.messages = suppress.messages)

  nca_df <- nca_df |> dplyr::select(-param_sd_z, -param_units) |>
    tidyr::pivot_wider(names_from = param_name,
                       values_from = param_value)

  if (win %in% TRUE) {
    nca_df <- nca_df |>
      dplyr::right_join(winmodel_df,
                        relationship = "many-to-many",
                        by = c(data_grp_vars)) |>
      dplyr::relocate(model, method, .after = Media)
  }


  # get tkstats
  tkstats_df <- get_tkstats.pk(obj = obj,
                               newdata = newdata,
                               model = model,
                               method = method,
                               tk_group = tk_group,
                               dose_norm = dose_norm,
                               exclude = exclude,
                               suppress.messages = suppress.messages) |>
    ungroup()

  if (win %in% TRUE) {
    tkstats_df <- dplyr::left_join(dplyr::ungroup(winmodel_df),
                                   tkstats_df,
                                   by = c(data_grp_vars, "method", "model"))
  }

  # prepare for merge
  nca_df_red <- nca_df |>
    dplyr::rename_with(~ paste0(.x, ".nca", recycle0 = TRUE),
                       !dplyr::any_of(c(grp_vars, "dose_norm", "model", "method")))

  tkstats_df_red <- tkstats_df |>
    dplyr::rename_with(~ paste0(.x, ".tkstats", recycle0 = TRUE),
                       !dplyr::any_of(c(grp_vars, "DATA_GROUP_ID",
                                             "Dose.Units", "Conc.Units",
                                             "Time.Units", "Time_trans.Units",
                                             "Rblood2plasma", "model", "method")
                       )
    )

  if (dose_norm) {
    message("eval_tkstats.pk(): ",
            "TK statistics AND NCA have been ",
            "calculated based on dose-normalized value of 1mg/kg"
    )
  }

  # Merge the tkstats and the nca data.frames
    tk_eval <- dplyr::left_join(
      tkstats_df_red,
      nca_df_red,
      by = c(grp_vars, "model", "method")
      ) |>
      dplyr::relocate(DATA_GROUP_ID, .before = 1L)

  # Filter out infinite and NA values for AUC_infinity
  if (finite_only) {
    tk_eval <- tk_eval |>
      filter(is.finite(AUC_infinity.tkstats),
             is.finite(AUC_infinity.nca))
  }

  return(tk_eval)
}

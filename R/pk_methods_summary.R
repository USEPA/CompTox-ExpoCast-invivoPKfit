#' Print summary of a `pk` object
#'
#' This summary includes summary information about the data; about any data
#' transformations applied; about the models being fitted; about the error model
#' being applied; and any fitting results, if the `pk` object has been fitted.
#' It also includes TK quantities calculated from the fitted model parameters,
#' e.g. halflife; clearance; tmax; Cmax; AUC; Css.
#'
#' @param object A [pk] object.
#' @param ... Additional arguments. Currently not in use.
#' @return A list of `data.frame`s consisting of a summary table of fitting options and results.
#' @export
#' @author Caroline Ring
#'
summary.pk <- function(object, ...) {
  objname <- deparse(substitute(object))
  status <- get_status(object)

  cat(paste0(objname, ":", "\n"))
  cat(paste0(attr(status, "msg"), "\n"))

  coef <- model <- method <- RMSE <- Fold_Error <- NULL
  avg_FoldErr <- median_FoldErr <- within_2fold <- NULL

  if (status < 5) {
    return(message("Please fit the object above before running, or try get_data_summary()"))
  }
    # things that require status 1
    # things that are just instructions:
    # data
    # preprocessing settings
    settings_preprocess <- get_settings_preprocess(object)
    cat(paste0(
      "\nData preprocessing settings:\n",
      paste(names(settings_preprocess),
            settings_preprocess,
            sep = " = ",
            collapse = "\n"),
      "\n"))

    # data info settings
    settings_data_info <- get_settings_data_info(object)
    cat(paste0(
      "\nData info settings:\n",
      paste(names(settings_data_info),
            settings_data_info,
            sep = " = ",
            collapse = "\n"),
      "\n"))
    # stat_model

    models <- names(object$stat_model)

    cat(paste0("\nModels to be fitted:\n",
               toString(models),
               "\n"))

    # stat_error_model

    error_group <- get_error_group(object)

    cat(paste0("\nError model grouping: ",
               rlang::as_label(error_group),
               "\n"))

    # optimx settings

    settings_optimx <- get_settings_optimx(object)
    cat(paste0("\nSettings for optimx::optimx():\n",
               paste(names(settings_optimx),
                     settings_optimx,
                     sep = " = ",
                     collapse = "\n"),
               "\n"))

  # scale_conc

  scale_conc <- get_scale_conc(object)
  cat(paste0(
    "\nConcentration scaling/transformation:\n",
    paste(names(scale_conc),
          scale_conc,
          sep = " = ",
          collapse = "\n"),
    "\n"))

  # scale_time

  scale_time <- get_scale_time(object)
  cat(paste0(
    "\nTime scaling/transformation:\n",
    paste(names(scale_time),
          scale_time,
          sep = " = ",
          collapse = "\n"),
    "\n"))

  # Produce a data.frame with all of the instructions in one place
  settings_preprocess_DF <- data.frame(instructions_category = "settings_preprocess",
                                       instruction_name = names(settings_preprocess),
                                       instruction_value = sapply(settings_preprocess,
                                                                  rlang::as_label))

  settings_data_info_DF <- data.frame(instructions_category = "settings_data_info",
                                      instruction_name = names(settings_data_info),
                                      instruction_value = sapply(settings_data_info,
                                                                 rlang::as_label))

  settings_optimx_DF <- data.frame(instructions_category = "settings_optimx",
                                   instruction_name = names(settings_optimx),
                                   instruction_value = sapply(settings_optimx,
                                                              rlang::as_label))

  stat_error_model_DF <- data.frame(instructions_category = "stat_error_model",
                                    instruction_name = "error_group",
                                    instruction_value = rlang::as_label(error_group))

  scale_conc_DF <- data.frame(instructions_category = "scale_conc",
                              instruction_name = names(scale_conc),
                              instruction_value = sapply(scale_conc,
                                                         rlang::as_label))

  scale_time_DF <- data.frame(instructions_category = "scale_time",
                              instruction_name = names(scale_time),
                              instruction_value = sapply(scale_time,
                                                         rlang::expr_deparse))

  instructions_DF <- do.call(rbind,
                             list(settings_preprocess_DF,
                                  settings_data_info_DF,
                                  scale_conc_DF,
                                  scale_time_DF,
                                  stat_error_model_DF,
                                  settings_optimx_DF))

  rownames(instructions_DF) <- NULL

  # things that require preprocessing & data_info (status 3):
  # data_summary
  cat("Data summary:\n")
  data_info <- get_data_info(object)
  print(data_info)
  cat("\n")

  # nca
  cat("NCA results:\n")
  nca <- get_nca(object)
  print(nca)
  cat("\n")

  # things that require pre-fitting (status 4):
  # stat_error_model
  message("Error Model Variables Specified: ")
  print(as.character(get_error_group(object))[-1])
  # stat_model par_DFs
  prefit_DF <- object$prefit$par_DF

  print(get_prefit(object))

  # things that require fitting (status 5):
  # coefficients
  # get model coefficients
  outDF <- suppressMessages(left_join(coef(object, include_NAs = TRUE),
                     coef_sd(object, table_format = TRUE)))

  # goodness of fit metrics
  aic_all <- suppressMessages(AIC(object)) # Already includes logLik
  bic_all <- suppressMessages(BIC(object))
  rmse_all <- rmse(object, use_scale_conc = FALSE) |>
    dplyr::group_by(!!!object$data_group, model, method, Route, Media, Dose) |>
    dplyr::mutate(avg_rmse = mean(RMSE, na.rm = TRUE)) |>
    tidyr::nest(full_rmse = c(Time, RMSE)) |>
    suppressMessages()

  fold_errors_all <- fold_error(object) |>
    dplyr::group_by(!!!object$data_group, model, method, Route, Media, Dose) |>
    dplyr::mutate(avg_FoldErr = mean(Fold_Error, na.rm = TRUE),
                  median_FoldErr = median(Fold_Error, na.rm = TRUE),
                  within_2fold = sum(Fold_Error >= 0.5 & Fold_Error <= 2) / length(Fold_Error)) |>
    dplyr::select(Time, Fold_Error,
                  avg_FoldErr, median_FoldErr, within_2fold) |>
    tidyr::nest(full_fold_errors = c(Time, Fold_Error)) |>
    suppressMessages()

  gof_DF <- dplyr::inner_join(AIC(object), BIC(object)) |>
    dplyr::right_join(dplyr::inner_join(rmse_all, fold_errors_all)) |>
    suppressMessages()


  # model comparison
  winmodel_DF <- suppressMessages(get_winning_model(object))

  # get TK stats
  tkstats <- suppressMessages(get_tkstats(object))
  # pivot wider and rowbind?


  # compare Cmax, AUC of winning model (for each method) to NCA Cmax, AUC, for each NCA group
  tk_nca <- suppressMessages(eval_tkstats(object))


  return(list("instructions" = instructions_DF,
              "data_summary" = get_data_info(object),
              "nca" = get_nca(object),
              "prefit_bounds" = prefit_DF,
              "parameters" = outDF,
              "goodness_of_fit" = gof_DF,
              "winning_model" = winmodel_DF,
              "tkstats" = tkstats,
              "tkstats_winning_vs_nca" = tk_nca))
}

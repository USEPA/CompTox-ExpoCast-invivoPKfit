#' Print a user-friendly version of a `pk` object
#'
#' Prints a clear summary of `pk` object and returns it invisibly.
#'
#' @param x A pk object.
#' @param ... Additional arguments. Currently not in use.
#'
#' @return Summary output
#' @export
#' @author Gilberto Padilla Mercado
print.pk <- function(x, ...) {

  status <- get_status(x, suppress.messages = FALSE)
  data_group_vars <- sapply(x$data_group, rlang::as_label)
  error_group_vars <- sapply(x$stat_error_model$error_group, rlang::as_label)
  models_present <- names(x$stat_model)
  scale_conc <- get_scale_conc(x)
  scale_time <- get_scale_time(x)
  summary_group_vars <- sapply(x$settings_data_info$summary_group, rlang::as_label)

  if (status > 1) {
    cat("Input data size:\n",
        "Rows: ", dim(x$data)[1], "\tColumns: ", dim(x$data)[2], "\n")
    cat("Unique Chemicals: ", length(unique(x$data$Chemical)), "\n")
    cat("Unique Species: ", length(unique(x$data$Species)), "\n")
    cat("Unique Routes: ", length(unique(x$data$Route)), "\n")
    cat("Unique Mediums: ", length(unique(x$data$Media)), "\n\n")
    if (all(error_group_vars %in% data_group_vars)) {
      cat("The error model is pooled.\n")
    } else {
      cat("The error model is hierarchical.\n")
    }
    cat("Data Group:\t", toString(data_group_vars), "\n")
    cat("Error Group:\t", toString(error_group_vars), "\n")
    cat("Scaling:\n\t",
        "Dose-normalization? ", scale_conc$dose_norm, "\n\t",
        "Log10-transformation? ", scale_conc$log10_trans, "\n\t")
    if (scale_time$new_units %in% "identity") {
      cat("Time scales are not adjusted.\n")
    } else if (scale_time$new_units %in% "auto") {
      cat("Time scales automatically adjusted.\n")
    } else {
      cat("Times converted to:\t", scale_time$new_units, "\n")
    }
    cat("Models: ", models_present, "\n\n")
  }

  if (status > 3) {
    fit_check_out <- x$prefit$fit_check
    cat("Average number of detected observations: ",
        round(mean(fit_check_out$n_detect, na.rm = TRUE) / length(models_present),
              digits = 0),
        "\n\n")
    cat("Fit decision summary:\n")
    cat("Continue:\t", sum(fit_check_out$fit_decision %in% "continue"), "\n")
    cat("Abort:\t\t", sum(fit_check_out$fit_decision %in% "abort"), "\n\n")
  }

  if (status == 5) {
    fit_results <- x$fit[c(data_group_vars,
                           "model", "method", "xtime", "convcode", "message")]
    fit_results <- unique(fit_results)
    cat("Fitted object summary:\n")
    cat("Average time per Data Group: ",
        signif(mean(fit_results$xtime, na.rm = TRUE), 3), "\n")
    cat("Number of non-zero convergence codes: ",
        sum(fit_results$convcode != 0), "\n\n")

  }

  invisible(x)
}

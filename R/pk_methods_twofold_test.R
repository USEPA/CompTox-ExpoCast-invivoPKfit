#' Evaluate whether data and predictions are within two-fold of mean or concentration, respectively
#'
#' At each timepoint across CvT experimental data, there are three ways that data may
#' be presented. These can be found as either:
#' - multiple individual observations
#' - single individual observation
#' - summarized group of observations (mean concentration and standard deviation)
#'
#' For the purposes of this calculations we largely divide the data into two groups,
#' those with individual observations, where N_Subjects == 1, and the summarized
#' group of observations.
#'
#' First this creates mean-normalized concentrations for individual data.
#' Then it summarizes data (individual & summarized) by `mean` and `sd`.
#' It tests whether predictions are within two-fold of mean,
#' in the latter case whether the 95% of observations are within 2-fold.
#'
#' Furthermore if `pk` object `status == 5` then it calculates the model error by
#' evaluating _prediction/concentration_ at each timepoint for all data. Each test
#' is done for data from individual subject observations and for all data by summarizing the
#' observations.
#'
#' Only non-excluded detects are included in this analysis.
#'
#' @param obj A pk object.
#' @param sub_pLOQ TRUE (default): Substitute all predictions below the LOQ with
#'   the LOQ before computing fold errors. FALSE: do not. Only used if `obj` has been fitted and predictions are possible.
#' @param suppress.messages Logical: whether to suppress message printing. If
#'   NULL (default), uses the setting in `obj$settings_preprocess$suppress.messages`.
#' @param model Optional: Specify one or more of the fitted models for which to
#'   make predictions. If NULL (the default), predictions will be returned for
#'   all of the models in `object$stat_model`.
#' @param method Optional: Specify one or more of the [optimx::optimx()] methods
#'   for which to make predictions. If NULL (the default), predictions will be
#'   returned for all of the models in `object$settings_optimx$method`.
#' @param ... Additional arguments. Currently unused.
#'
#' @return A list of data frames.
#' @export
#' @author Gilberto Padilla Mercado
#'
twofold_test.pk <- function(obj,
                            sub_pLOQ = TRUE,
                            suppress.messages = NULL,
                            model = NULL,
                            method = NULL,
                            ...) {

  if (is.null(suppress.messages)) {
    suppress.messages <- obj$settings_preprocess$suppress.messages
  }

  # Check the status of the pk object
  status <- suppressMessages(get_status(obj))
  if (status < 2) {
    stop("twofold_test.pk(): Please do data preprocessing step on this pk object",
         " by running do_preprocess(...)"
         )
  }

  # Get the data to evaluate
  data_cvt <- get_data(obj)

  # Only get observations where 'Detect' is TRUE and 'exclude' is FALSE
  # because other observations won't have credible variability from mean
  data_cvt <- subset(data_cvt, subset = (Detect | !exclude | Conc > LOQ))




  # Initialize output list
  out_list <- list()

  # Set variable describing necessary columns
  vital_col <- c(
    "Chemical",
    "Species",
    "Reference",
    "Route",
    "Media",
    "Dose",
    "Time",
    "Time.Units",
    "Conc",
    "Conc_SD",
    "Conc.Units",
    "N_Subjects",
    "Detect",
    "exclude"
  )

  if (!all(vital_col %in% names(data_cvt))) {
    stop("Missing the following columns:",
         toString(setdiff(vital_col, names(data_cvt)))
         )
  }

  data_cvt <- data_cvt[vital_col]


  # Find how many detected observations per Chemical + Species + Route + Media + Dose
  data_counts <- data_cvt |>
    dplyr::filter(Detect) |>
    dplyr::count(Chemical, Species, Reference,
                 Route, Media, Dose, Time, N_Subjects,
                 name = "Count")

  # There must be at least 2 values per timepoint for individual data or else
  # Conc is equal to mean(Conc)
  # Each Summarized Group of Observations should have multiple subjects
  # For individual data, can only summarize if there are multiple observations per timepoint
  # Note we only use data_counts filtering left join with individual data
  # Group data have multiple subjects per experimental timepoint
  data_counts <- dplyr::left_join(data_counts, data_cvt,
                       by = c("Chemical", "Species", "Reference",
                              "Route", "Media", "Dose", "Time",
                              "N_Subjects")) |>
    dplyr::mutate(data_descr = dplyr::case_when(
      (Count > 1 & N_Subjects == 1) ~ "Individual Data, Multiple Observations",
      (Count == 1 & N_Subjects == 1) ~ "Individual Data, Single Observation",
      (Count == 1 & N_Subjects > 1 & Conc_SD > 0) ~ "Grouped Data w SD",
      (N_Subjects > 1 & Conc_SD == 0) ~ "Grouped Data no SD",
      .default = NA_character_
    ))

  out_list[["data_descriptors"]] <- data_counts

  # There must be at least 2 values per timepoint for individual data or else
  # Conc is equal to mean(Conc)
  indiv_data <- subset(data_counts, subset = (data_descr %in% "Individual Data, Multiple Observations"))
  # I will also save the individual single observations
  single_data <- subset(data_counts, subset = (data_descr %in% "Individual Data, Single Observation"))
  # Grouped data
  sgroup_data <- subset(data_counts, subset = (data_descr %in% "Grouped Data w SD"))
  sgroup_nsd_data <- subset(data_counts, subset = (data_descr %in% "Grouped Data no SD"))

  sgroup_data <- dplyr::mutate(sgroup_data,
                               conc_mean = Conc, conc_sd = Conc_SD)

  sgroup_nsd_data <- dplyr::mutate(sgroup_nsd_data,
                                   conc_mean = Conc, conc_sd = Conc_SD)

  # Next I will deal with individual data and
  # calculate mean and standard deviation for individual data
  indiv_data_summary <- indiv_data |>
    dplyr::group_by(Chemical, Species, Reference, Route, Media, Dose, Time) |>
    dplyr::summarise(conc_mean = mean(Conc, na.rm = TRUE),
                     conc_sd = sd(Conc, na.rm = TRUE))

  # Need to ensure both data.frames have the same colnames before dplyr::bind_rows
  indiv_data_summary <- indiv_data_summary[c(intersect(names(indiv_data_summary),
                                                       names(sgroup_data)))]
  sgroup_data <- sgroup_data[c(intersect(names(indiv_data_summary),
                                       names(sgroup_data)))]

  # Combined summarized data (indiv + sgroup)
  total_data_summary <- dplyr::bind_rows(indiv_data_summary, sgroup_data) |>
    dplyr::ungroup() |>
    dplyr::mutate(twofold_95 = ((conc_mean + 2 * conc_sd) / conc_mean) <= 2)

  twofold_95 <- total_data_summary |>
    dplyr::group_by(Route) |>
    dplyr::summarise(within = sum(twofold_95),
                     outside = sum(!twofold_95))
  all_95 <- total_data_summary |>
    dplyr::summarise(within = sum(twofold_95),
                     outside = sum(!twofold_95)) |>
    dplyr::mutate(Route = "All", .before = 1)

  twofold_95 <- dplyr::bind_rows(twofold_95, all_95)

  ### General Test Summary of Data that are NOT single observation per timepoint
  twofold_95 <- rowwise_calc_percentages(twofold_95, group_cols = "Route")
  ###

  # Use indiv_data (id_) for individual data evaluations
  # Use total_data_summary (tds_) for summarized data evaluations
  # Add grouped mean of individual data to indiv_data
  id_mean <- suppressMessages(dplyr::left_join(indiv_data, indiv_data_summary)) |>
    dplyr::mutate(foldConc = Conc / conc_mean)

  # Make a table of values with percent within or outside factors of two
  # Summarizing number of fold concentrations from the mean
  id_twofold_route <- id_mean |>
    dplyr::group_by(Route) |>
    dplyr::summarise(above_twofold = sum(foldConc > 2),
                     within_twofold = sum(foldConc >= 0.5 & foldConc <= 2),
                     below_twofold = sum(foldConc < 0.5))

  # Also need a row for total
  id_total_twofold <- id_mean |>
    dplyr::summarise(above_twofold = sum(foldConc > 2),
                     within_twofold = sum(foldConc >= 0.5 & foldConc <= 2),
                     below_twofold = sum(foldConc < 0.5)) |>
    dplyr::mutate(Route = "All", .before = 1)

  ### Combine into single data.frame
  id_total_twofold <- dplyr::bind_rows(id_twofold_route, id_total_twofold)

  id_total_twofold <- rowwise_calc_percentages(id_total_twofold,
                                               group_cols = "Route")


  # First three outputs
  out_list[["individual_data"]] <- id_mean
  out_list[["summarized_data_test"]] <- twofold_95
  out_list[["indiv_data_test"]] <- id_total_twofold

# If predictions are possible
  if (status == 5) {

    winmodel <- get_winning_model(obj = obj)

    # Only keep detects and non-excluded (observations to compare to)
    pred_win <- fold_error(obj = obj,
                 model = model,
                 method = method,
                 sub_pLOQ = sub_pLOQ,
                 suppress.messages = suppress.messages) |>
        dplyr::ungroup() |>
        dplyr::filter(Detect, !exclude) |>
        dplyr::semi_join(winmodel) |>
      suppressMessages()

    # Make a table of values with percent within or outside factors of two
    # Summarizing number of fold concentrations from the mean
    # Here I really only care about 1- and 2-compartment models
    pred_twofold_summary <- pred_win |>
      dplyr::group_by(model, method) |>
      dplyr::summarise(above_twofold = sum(Fold_Error > 2),
                       within_twofold = sum(dplyr::between(Fold_Error, 0.5, 2)),
                       below_twofold = sum(Fold_Error < 0.5))

    ###
    # Get totals for method
    all_models <- pred_win |>
      dplyr::filter(!(model %in% "model_flat")) |>
      dplyr::group_by(method) |>
      dplyr::summarise(above_twofold = sum(Fold_Error > 2),
                       within_twofold = sum(dplyr::between(Fold_Error, 0.5, 2)),
                       below_twofold = sum(Fold_Error < 0.5)) |>
      dplyr::mutate(model = "All non-flat", .before = 1)

    full_pred_twofold <- dplyr::bind_rows(pred_twofold_summary,
                               all_models)

    full_pred_twofold <- rowwise_calc_percentages(full_pred_twofold,
                                                  group_cols = c("model",
                                                                 "method"))

    ##
    # Output these data.frames
    out_list[["model_error_all"]] <- pred_win
    out_list[["model_error_summary"]] <- full_pred_twofold

    ##
    # This merge should have both Fold_Error AND foldConc
    #inner join keeps only observations in both tables,
    #keyed by all variables in common between pred_win and id_mean,
    #so it's only keeping the time points with replicates
    data_preds <- suppressMessages(
      dplyr::inner_join(pred_win,
                        id_mean,
                        relationship = "many-to-many")
    )

    data_preds["both_within"] <- (
      dplyr::between(data_preds$Fold_Error, 0.5, 2) &
        dplyr::between(data_preds$foldConc, 0.5, 2)
    )

    data_preds["both_outside"] <- (
      !dplyr::between(data_preds$Fold_Error, 0.5, 2) &
        !dplyr::between(data_preds$foldConc, 0.5, 2)
    )

    data_preds["model_outside"] <- (
      !dplyr::between(data_preds$Fold_Error, 0.5, 2) &
        dplyr::between(data_preds$foldConc, 0.5, 2)
    )

    data_preds["data_outside"] <- (
      dplyr::between(data_preds$Fold_Error, 0.5, 2) &
        !dplyr::between(data_preds$foldConc, 0.5, 2)
    )
    # Make a table of values with percent within or outside factors of two
    # Summarizing number of fold concentrations from the mean
    data_pred_summary <- data_preds |>
      dplyr::group_by(model, method) |>
      dplyr::summarise(across(c(.data$both_within:.data$data_outside), sum))
    data_pred_summary_combo <- data_preds |>
      dplyr::group_by(method) |>
      dplyr::summarise(across(c(.data$both_within:.data$data_outside), sum)) |>
      dplyr::mutate(model = "All", .before = 1)
    data_pred_summary <- dplyr::bind_rows(data_pred_summary, data_pred_summary_combo)

    data_pred_summary <- rowwise_calc_percentages(data_pred_summary,
                                                  group_cols = c("model",
                                                                 "method"))

    out_list[["indiv_data_fold_errors"]] <- data_preds
    out_list[["indiv_data_test_fold_errors"]] <- data_pred_summary
    out_list[["fold_errors_data_descriptors"]] <- merge(pred_win, data_counts,
                                                          all.x = TRUE)

    ###
    # Here is when we take the 95% within two-fold summary data into account
    data95_preds <- dplyr::left_join(pred_win,
                                     total_data_summary) |>
      dplyr::filter(!is.na(twofold_95))

    data95_preds["both_within"] <- (dplyr::between(data95_preds$Fold_Error, 0.5, 2) &
                                      (data95_preds$twofold_95))

    data95_preds["both_outside"] <- (!dplyr::between(data95_preds$Fold_Error, 0.5, 2) &
                                       !(data95_preds$twofold_95))

    data95_preds["model_outside"] <- (!dplyr::between(data95_preds$Fold_Error, 0.5, 2) &
                                        (data95_preds$twofold_95))
    data95_preds["data_outside"] <- (dplyr::between(data95_preds$Fold_Error, 0.5, 2) &
                                       !(data95_preds$twofold_95))

    # Make a table of values with percent within or outside factors of two
    # Summarizing number of fold concentrations from the mean
    data95_pred_summary <- data95_preds |>
      dplyr::group_by(model, method) |>
      dplyr::summarise(across(c(.data$both_within:.data$data_outside), sum))

    data95_pred_summary <- rowwise_calc_percentages(data95_pred_summary,
                                                  group_cols = c("model",
                                                                 "method"))

    out_list[["summarized_data_model_error"]] <- data95_pred_summary

    ###
  } else { # This is the case where the model has not been fit
    ###
    warning("pk object not fit, unable to run metrics for predictions")
  }

  if (!suppress.messages) {
    cat("\n\n")
    message("Data Summary of 95% of mean-normalized Concentrations within 2-fold Test")
    print(out_list[["summarized_data_test"]])

    # Print the summaries
    if (suppress.messages %in% FALSE) {
    cat("\n\n")
    message("Result: Individual Data Fold Concentration from Mean 2-fold Test")
    print(out_list[["indiv_data_test"]])

    if (status == 5) {
      cat("\n\n")
      message("Result: Model Fold Error Summary")
      print(out_list[["model_error_summary"]])

      cat("\n\n")
      message("Result: Individual Data and Model Error Test")
      print(out_list[["indiv_data_test_fold_errors"]])
    }
    }
  }

  invisible(out_list)

}


#' Helper function for calculating percentages of count data, by row
#'
#' This function takes totals and calculates rowise percentages across columns
#' Expects columns for each percentage, can specify a vector of "grouping" column names
#'
#' @param data A data.frame that contains columns of count data and
#' possibly columns of group names.
#' @param group_cols String or numeric indices for the columns which contain
#' grouping variables.
#'
#' @return data.frame with rowwise totals and percentages
rowwise_calc_percentages <- function(data,
                                     group_cols = NULL) {
  if (!is.null(group_cols)) {
    if (is.list(group_cols)) group_cols <- unlist(group_cols)
    group_cols <- unique(group_cols)
    op_cols <- setdiff(names(data), names(data[group_cols]))
  } else {
    op_cols <- names(data)
  }

  # Calculate total
  data$total <- rowSums(data[op_cols])

  # Calculate percentages
  data[paste0("percent_", op_cols)] <- apply(
    data[op_cols], 2,
    \(x) signif(100 * x / data$total, 4)
  )

  return(data)
}

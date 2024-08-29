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
#' @param obj A pk object.
#' @param suppress_messages Logical, default: FALSE.
#' Should there be no printed output?
#'
#' @return A list of data frames.
#' @export
#' @author Gilberto Padilla Mercado
#'
twofold_test.pk <- function(obj,
                            suppress_messages = FALSE) {
  # Check the status of the pk object
  status <- suppressMessages(get_status(obj))
  if (status < 2) {
    stop(paste0("Please do data preprocessing step on this pk object",
                " by running do_preprocess(...)"))
  }

  # Get the data to evaluate
  data_cvt <- get_data(obj)

  # Only get observations where 'Detect' is TRUE and 'exclude' is FALSE
  # because these observations won't have credible variability from mean
  data_cvt <- subset(data_cvt,
                    subset = (Detect | !exclude | Conc > LOQ))


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
    stop(paste("Missing the following columns:",
               paste(setdiff(vital_col, names(data_cvt)),
                     collapse = ", ")))
  }

  data_cvt <- data_cvt[vital_col]


  # Find how many observations per Chemical + Species + Route + Media + Dose
  data_counts <- aggregate(
    Conc ~ Chemical + Species + Reference + Route + Media + Dose + Time + N_Subjects,
    data = data_cvt,
    FUN = NROW)
  names(data_counts)[names(data_counts) == 'Conc'] <- 'Count'

  # There must be at least 2 values per timepoint for individual data or else
  # Conc = mean(Conc)
  # Each Summarized Group of Observations should have multiple subjects
  # For individual data, can only summarize if there are multiple observations per timepoint
  # Note we only use data_counts filtering left join with individual data
  # Group data have multiple subjects per experimental timepoint
  data_counts <- merge(data_counts, data_cvt, all.x = TRUE,
                       by = c("Chemical", "Species", "Reference",
                              "Route", "Media", "Dose", "Time",
                              "N_Subjects"))
  data_counts$data_descr <- with(data_counts,
                                 dplyr::case_when(
                                   (Count > 1 & N_Subjects == 1) ~ "Individual Data Multiple Obs",
                                   (Count == 1 & N_Subjects == 1) ~ "Individual Data Single Obs",
                                   (Count == 1 & N_Subjects > 1 & Conc_SD > 0) ~ "Grouped Data w SD",
                                   (N_Subjects > 1 & Conc_SD == 0) ~ "Grouped Data no SD",
                                   .default = "Uncharacterized"
  ))

  out_list[['data_descriptors']] <- data_counts

  # There must be at least 2 values per timepoint for individual data or else
  # Conc = mean(Conc)
  indiv_data <- subset(data_counts, subset = (data_descr %in% "Individual Data Multiple Obs"))
  # I will also save the individual single obsevations
  single_data <- subset(data_counts, subset = (data_descr %in% "Individual Data Single Obs"))
  # Grouped data
  sgroup_data <- subset(data_counts, subset = (data_descr %in% "Grouped Data w SD"))
  sgroup_nsd_data <- subset(data_counts, subset = (data_descr %in% "Grouped Data noSD"))


  sgroup_data['conc_mean'] <- sgroup_data['Conc']
  sgroup_data['conc_sd'] <- sgroup_data['Conc_SD']

  sgroup_nsd_data['conc_mean'] <- sgroup_nsd_data['Conc']
  sgroup_nsd_data['conc_sd'] <- sgroup_nsd_data['Conc_SD']


  # Next I will deal with individual data and
  # calculate mean and standard deviation for individual data
  indiv_data_summary <- do.call(data.frame,
                        aggregate(
                          Conc ~ Chemical + Species + Reference + Route + Media + Dose + Time,
                          data = indiv_data,
                          FUN = \(x) {
                            c(
                              conc_mean = mean(x, na.rm = TRUE),
                              conc_sd = sd(x, na.rm = TRUE)
                            )
                          }))
  # Rename new columns (aggregate adds "Var." for the "Var ~ " column being summarized)
  names(indiv_data_summary) <- gsub(pattern = "Conc.",
                            replacement = "",
                            x = names(indiv_data_summary))

  # Need to ensure both data.frames have the same colnames before rbind
  indiv_data_summary <- indiv_data_summary[c(intersect(names(indiv_data_summary),
                                                       names(sgroup_data)))]
  sgroup_data <- sgroup_data[c(intersect(names(indiv_data_summary),
                                       names(sgroup_data)))]

  # Combined summarized data (indiv + sgroup)
  total_data_summary <- rbind(indiv_data_summary, sgroup_data)

  total_data_summary['twofold_95'] <- with(total_data_summary,
                                   ifelse((conc_mean + (2*conc_sd))/conc_mean <= 2,
                                          TRUE, FALSE))

  twofold_95 <- do.call(data.frame,
                        aggregate(twofold_95 ~ Route,
                                  data = total_data_summary,
                                  FUN = \(x) {
                                    c(within = sum(x),
                                      outside = sum(!x))
                                  }))
  names(twofold_95) <- gsub(pattern = "\\.",
                                  replacement = "_",
                                  x = names(twofold_95))

  twofold_95 <- rbind(twofold_95,
                      data.frame(Route = "All",
                                 twofold_95_within = with(twofold_95,
                                                          sum(twofold_95_within)),
                                 twofold_95_outside = with(twofold_95,
                                                           sum(twofold_95_outside))
                      )
  )
  ### General Test Summary of Data that are NOT single observation per timepoint
  twofold_95 <- rowwise_calc_percentages(twofold_95, group_cols = "Route")
  ###

  # Use indiv_data (id_) for individual data evaluations
  # Use total_data_summary (tds_) for summarized data evaluations
  # Add grouped mean of individual data to indiv_data
  id_mean <- merge(indiv_data, indiv_data_summary,
                     all.x = TRUE)

  # Calculate the fold concentration from the mean
  id_mean['foldConc'] <- with(id_mean, Conc / conc_mean)

  # Make a table of values with percent within or outside factors of two
  # Summarizing number of fold concentrations from the mean
  id_twofold_route <- do.call(data.frame,
                                  aggregate(foldConc ~ Route,
                                            data = id_mean,
                                            FUN = \(x) {
                                              c(
                                                above_twofold = sum(x > 2),
                                                within_twofold = sum(x >= 0.5 & x <= 2),
                                                below_twofold = sum(x < 0.5)
                                              )
                                            }
                                  ))
  names(id_twofold_route) <- gsub(pattern = "foldConc.",
                                    replacement = "",
                                    x = names(id_twofold_route))


  # Also need a row for total
  id_total_twofold <- data.frame(
    Route = "All",
    above_twofold = with(id_twofold_route,
                         sum(above_twofold)),
    within_twofold = with(id_twofold_route,
                          sum(within_twofold)),
    below_twofold = with(id_twofold_route,
                         sum(below_twofold))
  )

  ### Combine into single data.frame
  id_total_twofold <- rbind(id_twofold_route, id_total_twofold)

  id_total_twofold <- rowwise_calc_percentages(id_total_twofold,
                                               group_cols = "Route")


  # First two outputs
  out_list[["summarized_data_test"]] <- twofold_95
  out_list[['indiv_data_test']] <- id_total_twofold

# If predictions are possible
  if (status == 5) {
    winmodel <- get_winning_model(obj = obj)

    pred_win <- merge(winmodel,
                      suppressMessages(fold_errors(obj = obj,
                                                   type = "conc"))
    )

    # Make a table of values with percent within or outside factors of two
    # Summarizing number of fold concentrations from the mean
    # Here I really only care about 1- and 2-compartment models
    pred_twofold_summary <- do.call(
      data.frame,
      aggregate(Fold_Error ~ Route + model + method,
                data = subset(pred_win,
                              subset = (model %in% c("model_1comp", "model_2comp"))),
                FUN = \(x) {
                  c(
                    above_twofold = sum(x > 2),
                    within_twofold = sum(x >= 0.5 & x <= 2),
                    below_twofold = sum(x < 0.5)
                  )
                }
      ))

    names(pred_twofold_summary) <- gsub(pattern = "Fold_Error.",
                                        replacement = "",
                                        x = names(pred_twofold_summary))

    ###
    # Get totals for models + method
    all_routes <- cbind(
      Route = "All",
      merge(
        merge(
          aggregate(above_twofold ~ model + method,
                    pred_twofold_summary, sum),
          aggregate(within_twofold ~ model + method,
                    pred_twofold_summary, sum)
        ),
        aggregate(below_twofold ~ model + method,
                  pred_twofold_summary, sum)
      ))


    ###
    # Get totals for method
    all_models <- cbind(
      Route = "All",
      model = "1- & 2-comp",
      merge(
        merge(
          aggregate(above_twofold ~ method,
                    pred_twofold_summary, sum),
          aggregate(within_twofold ~ method,
                    pred_twofold_summary, sum)
        ),
        aggregate(below_twofold ~ method,
                  pred_twofold_summary, sum)
      ))

    full_pred_twofold <- rbind(pred_twofold_summary,
                               all_routes,
                               all_models)

    full_pred_twofold <- rowwise_calc_percentages(full_pred_twofold,
                                                  group_cols = c("Route",
                                                                 "model",
                                                                 "method"))

    ##
    # Output these data.frames
    out_list[['model_error_twofold_summary']] <- full_pred_twofold

    ##
    # This merge should have both Fold_Error AND foldConc
    data_preds <- merge(id_mean, pred_win, all.x = TRUE)

    data_preds['both_within'] <- (dplyr::between(data_preds$Fold_Error, 0.5, 2) &
                                    dplyr::between(data_preds$foldConc, 0.5, 2))

    data_preds['both_outside'] <- (!dplyr::between(data_preds$Fold_Error, 0.5, 2) &
                                     !dplyr::between(data_preds$foldConc, 0.5, 2))

    data_preds['model_outside'] <- (!dplyr::between(data_preds$Fold_Error, 0.5, 2) &
                                      dplyr::between(data_preds$foldConc, 0.5, 2))
    data_preds['data_outside'] <- (dplyr::between(data_preds$Fold_Error, 0.5, 2) &
                                     !dplyr::between(data_preds$foldConc, 0.5, 2))

    # Make a table of values with percent within or outside factors of two
    # Summarizing number of fold concentrations from the mean
    data_pred_summary <- do.call(data.frame,
                                 aggregate(
                                   cbind(both_within,
                                         both_outside,
                                         model_outside,
                                         data_outside) ~ model + method,
                                   data = data_preds,
                                   FUN = \(x) {sum(x + 0)}))

    data_pred_summary <- rbind(
      data_pred_summary,
      cbind(
        model = "1- & 2-comp",
        merge(
          merge(
            merge(
              aggregate(both_within ~ method,
                        subset(data_pred_summary,
                               subset = (model %in% c("model_1comp",
                                                      "model_2comp"))), sum),
              aggregate(both_outside ~ method,
                        subset(data_pred_summary,
                               subset = (model %in% c("model_1comp",
                                                      "model_2comp"))), sum)
            ),
            aggregate(model_outside ~ method,
                      subset(data_pred_summary,
                             subset = (model %in% c("model_1comp",
                                                    "model_2comp"))), sum)
          ),
          aggregate(data_outside ~ method,
                    subset(data_pred_summary,
                           subset = (model %in% c("model_1comp",
                                                  "model_2comp"))), sum)
        )
      )
    )


    data_pred_summary <- rowwise_calc_percentages(data_pred_summary,
                                                  group_cols = c("model",
                                                                 "method"))

    out_list[['indiv_data_fold_errors']] <- data_preds
    out_list[['indiv_data_test_fold_errors']] <- data_pred_summary
    out_list[["fold_errors_data_descriptors"]] <- merge(pred_win, data_counts,
                                                          all.x = TRUE)

    ###
    # Here is when we take the 95% within two-fold summary data into account
    data95_preds <- merge(total_data_summary, pred_win, all.x = TRUE)

    data95_preds['both_within'] <- (dplyr::between(data95_preds$Fold_Error, 0.5, 2) &
                                      (data95_preds$twofold_95))

    data95_preds['both_outside'] <- (!dplyr::between(data95_preds$Fold_Error, 0.5, 2) &
                                       !(data95_preds$twofold_95))

    data95_preds['model_outside'] <- (!dplyr::between(data95_preds$Fold_Error, 0.5, 2) &
                                        (data95_preds$twofold_95))
    data95_preds['data_outside'] <- (dplyr::between(data95_preds$Fold_Error, 0.5, 2) &
                                       !(data95_preds$twofold_95))

    # Make a table of values with percent within or outside factors of two
    # Summarizing number of fold concentrations from the mean
    data95_pred_summary <- do.call(data.frame,
                                 aggregate(
                                   cbind(both_within,
                                         both_outside,
                                         model_outside,
                                         data_outside) ~ model + method,
                                   data = data95_preds,
                                   FUN = \(x) {sum(x + 0)}))

    data95_pred_summary <- rowwise_calc_percentages(data95_pred_summary,
                                                  group_cols = c("model",
                                                                 "method"))

    out_list[['summarized_data_model_error']] <- data95_pred_summary

    ###
  } else { # This is the case where the model has not been fit
    ###
    warning("pk object not fit, unable to run metrics for predictions")
  }

  if (!suppress_messages) {
    cat("\n\n")
    message("Data Summary of 95% of mean-normalized Concentrations within 2-fold Test")
    print(out_list[['summarized_data_test']])

    # Print the summaries
    cat("\n\n")
    message("Result: Individual Data Fold Concentration from Mean 2-fold Test")
    print(out_list[['indiv_data_test']])

    if (status == 5) {
      cat("\n\n")
      message("Result: Model Fold Error Summary")
      print(out_list[['model_error_twofold_summary']])

      cat("\n\n")
      message("Result: Individual Data and Model Error Test")
      print(out_list[['indiv_data_test_fold_error']])

      cat("\n\n")
      message("Result: Summary Data and Model Error Test")
      print(out_list[['summarized_data_model_error']])
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

  data$total <-  apply(data[op_cols], 1, sum)

  # Calculate percentages
  data[paste0("percent_", op_cols)] <- apply(data[op_cols], 2,
                                             \(x) signif(100 * x / data$total,
                                                         4))

  return(data)
}

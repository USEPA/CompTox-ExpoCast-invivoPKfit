#' Evaluate whether data and predictions are within two-fold of mean or concentration, respectively
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
    Conc ~ Chemical + Species + Reference + Route + Media + Dose + Time,
    data = data_cvt,
    FUN = NROW)
  names(data_counts)[names(data_counts) == 'Conc'] <- 'Count'

  # There must be at least 2 values per timepoint for individual data
  data_counts <- subset(data_counts, subset = Count > 1)


  # For individual data, can only summarize if there are multiple observations per timepoint
  # Note we only use data_counts filtering left join with individual data
  indiv_data <- merge(data_counts[names(data_counts) != 'Count'],
                      subset(data_cvt, subset = ((N_Subjects == 1 | is.na(N_Subjects)) &
                                                  Conc_SD == 0)),
                      all.x = TRUE)

  # Group data have multiple subjects per experimental timepoint
  # Must have some data variability described
  group_data <- subset(data_cvt, subset = (N_Subjects >= 2 & Conc_SD > 0))
  group_data['conc_mean'] <- group_data['Conc']
  group_data['conc_sd'] <- group_data['Conc_SD']

  # Calculate mean and standard deviation for individual data
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
                                                       names(group_data)))]
  group_data <- group_data[c(intersect(names(indiv_data_summary),
                                       names(group_data)))]

  # Combined data
  total_data_summary <- rbind(indiv_data_summary, group_data)

  total_data_summary['twofold_95'] <- with(total_data_summary,
                                   ifelse((conc_mean + (2*conc_sd))/conc_mean <= 2,
                                          TRUE, FALSE))

  tds_list <- data.frame(
    Data = "Individual and Summarized Data",
    within_twofold_95 = 100*sum(total_data_summary['twofold_95'])/nrow(total_data_summary),
    outside_twofold_95 = 100*sum(!total_data_summary['twofold_95'])/nrow(total_data_summary))

  # Use indiv_data (id_) for individual data evaluations
  # Use total_data_summary (tds_) for summarized data evaluations


  # Add grouped mean of individual data to
  id_mean <- merge(indiv_data_summary, indiv_data,
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

  id_twofold_route <- twofold_calc_percentages(id_twofold_route)


  # Also need a row for total
  id_total_twofold <- data.frame(
    Route = "Both",
    above_twofold = with(id_twofold_route,
                         sum(above_twofold)),
    within_twofold = with(id_twofold_route,
                          sum(within_twofold)),
    below_twofold = with(id_twofold_route,
                         sum(below_twofold))
  )

  id_total_twofold <- twofold_calc_percentages(id_total_twofold)


  # First three outputs
  out_list[["Summarized Data Normal Test"]] <- tds_list
  out_list[['Individual Data by Route']] <- id_twofold_route
  out_list[['Individual Data All']] <- id_total_twofold



# If predictions are possible
  if (status == 5) {
    winmodel <- suppressMessages(get_winning_model(obj = obj))

    # id_mean already has foldConc to test
    # but needs exclude and Time.Units from initial data
    id_mean <- merge(id_mean,
                     data_cvt[vital_col])

    pred_twofold <- merge(winmodel,
                          predict(obj = obj,
                                  newdata = id_mean,
                                  use_scale_conc = FALSE,
                                  by_timepoint = FALSE,
                                  type = "conc"))


    pred_twofold['foldPred'] <- with(pred_twofold,
                                     Conc_est / Conc)



    # Make a table of values with percent within or outside factors of two
    # Summarizing number of fold concentrations from the mean
    pred_twofold_summary <- do.call(data.frame,
                                    aggregate(foldPred ~ Route + model + method,
                                              data = pred_twofold,
                                              FUN = \(x) {
                                                c(
                                                  above_twofold = sum(x > 2),
                                                  within_twofold = sum(x >= 0.5 & x <= 2),
                                                  below_twofold = sum(x < 0.5)
                                                )
                                              }
                                    ))

    names(pred_twofold_summary) <- gsub(pattern = "foldPred.",
                                        replacement = "",
                                        x = names(pred_twofold_summary))
    pred_twofold_summary <- twofold_calc_percentages(pred_twofold_summary)

    ###
    # Get totals for models + method
    both_routes <- cbind(
      merge(
        merge(
          aggregate(above_twofold ~ model + method,
                    pred_twofold_summary, sum),
          aggregate(within_twofold ~ model + method,
                    pred_twofold_summary, sum)
        ),
        aggregate(below_twofold ~ model + method,
                  pred_twofold_summary, sum),
      ))
    both_routes <- twofold_calc_percentages(both_routes)


    ###
    # Get totals for method
    all_models <- cbind(
      model = "All",
      merge(
        merge(
          aggregate(above_twofold ~ method,
                    pred_twofold_summary, sum),
          aggregate(within_twofold ~ method,
                    pred_twofold_summary, sum)
        ),
        aggregate(below_twofold ~ method,
                  pred_twofold_summary, sum),
      ))

    all_models <- twofold_calc_percentages(all_models)

    ##
    # Output these data.frames
    out_list[['Model Error by Route']] <- pred_twofold_summary
    out_list[['Model Error by Method and Model']] <- both_routes
    out_list[['Model Error by Method']] <- all_models

    browser()

    ##
    # This merge should have both foldPred AND foldConc
    # However, it also has multiple
    data_preds <- merge(pred_twofold,
                        data_mean)

    data_preds['both_within'] <- ((data_preds$foldPred >= 0.5 &
                                     data_preds$foldPred <= 2) &
                                    (data_preds$foldConc >= 0.5 &
                                       data_preds$foldConc <= 2))

    data_preds['both_outside'] <- ((data_preds$foldPred < 0.5 |
                                      data_preds$foldPred > 2) &
                                     (data_preds$foldConc < 0.5 |
                                        data_preds$foldConc > 2))

    data_preds['model_outside'] <- ((data_preds$foldPred < 0.5 |
                                       data_preds$foldPred > 2) &
                                      (data_preds$foldConc >= 0.5 |
                                         data_preds$foldConc <= 2))
    data_preds['data_outside'] <- ((data_preds$foldPred >= 0.5 |
                                      data_preds$foldPred <= 2) &
                                     (data_preds$foldConc < 0.5 |
                                        data_preds$foldConc > 2))

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

    data_pred_summary <- twofold_calc_percentages(data_pred_summary)


    ###
  } else { # This is the case where the model has not been fit
    ###
    warning("pk object not fit, unable to run metrics for predictions")
  }

  out_list[['Data with foldConc and foldPreds']] <- data_preds
  out_list[['foldPreds only']] <- pred_twofold

  if (!suppress_messages) {
    cat("\n\n")
    message("Data Summary of 95% of mean-normalized Concentrations within 2-fold")
    print(out_list[['data_95_twofold']])

    # Print the summaries
    cat("\n\n")
    message("Data Summary of Fold Concentration from Mean")
    print(out_list[['data_fold_summary']])

    if (status == 5) {
      cat("\n\n")
      message("Prediction / Concentration Summary")
      print(out_list[['pred_fold_summary']])
    }
  }

  return(out_list)
}



# This function takes totals and calculates rowise percentages across columns
twofold_calc_percentages <- function(data) {
  if (all(
    c("above_twofold",
      "within_twofold",
      "below_twofold") %in% names(data))) {

    data['total_twofold'] <-  with(data,
                                   above_twofold + within_twofold + below_twofold)
    # Calculate percentages
    data['percent_above'] <- with(data,
                                  signif(
                                    100 * above_twofold / total_twofold,
                                    4))
    data['percent_below'] <- with(data,
                                  signif(
                                    100 * below_twofold / total_twofold,
                                    4))
    data['percent_within'] <- with(data,
                                   signif(
                                     100 * within_twofold / total_twofold,
                                     4))

    return(data)
  } else if (all(
    c("both_within", "both_outside",
      "model_outside", "data_outside") %in% names(data)
  )) {
    data['total_twofold'] <-  with(data,
                                   both_within + both_outside + model_outside + data_outside)
    # Calculate percentages
    data['percent_both_win'] <- with(data,
                                     signif(
                                       100 * both_within / total_twofold,
                                       4))
    data['percent_both_out'] <- with(data,
                                     signif(
                                       100 * both_outside / total_twofold,
                                       4))
    data['percent_model_out'] <- with(data,
                                      signif(
                                        100 * model_outside / total_twofold,
                                        4))
    data['percent_data_out'] <- with(data,
                                     signif(
                                       100 * data_outside / total_twofold,
                                       4))
    return(data)
  } else {
    stop("The are some missing columns!")
  }
}


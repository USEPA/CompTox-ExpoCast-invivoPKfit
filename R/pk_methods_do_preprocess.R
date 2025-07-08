#' Do pre-processing
#'
#' Pre-process data for a `pk` object
#'
#' Data pre-processing for an object `obj` includes the following steps, in order:
#' \itemize{
#' \item Coerce data to class `data.frame` (if it is not already)
#' \item Rename variables to harmonized "`invivopkfit` aesthetic" variable names, using `obj$mapping`
#' \item Check that the data includes only routes in `obj$settings_preprocess$routes_keep` and media in `obj$settings_preprocess$media_keep`
#' \item Check that the data includes only one unit for concentration, one unit for time, and one unit for dose.
#' \item Coerce `Value`, `Value_SD`, `LOQ`, `Dose`, and `Time` to numeric, if they are not already.
#' \item Coerce `Species`, `Route`, and `Media` to lowercase.
#' \item Replace any negative `Value`, `Value_SD`, `Dose`, or `Time` with `NA`
#' \item If any non-NA `Value` is currently less than its non-NA LOQ, then replace it with NA
#' \item  Impute any NA `LOQ`: as `calc_loq_factor` * minimum non-NA `Value` in each `loq_group`
#' \item For any cases where `N_Subject`s is NA, impute `N_Subjects` = 1
#' \item For anything with `N_Subjects` == 1, set `Value_SD` to 0
#' \item Impute missing `Value_SD` as follows: For observations with `N_Subjects` > 1, take the minimum non-issing `Value_SD` for each `sd_group`. If all SDs are missing in an `sd_group`, then `Value_SD` for each observation in that group will be imputed as 0.
#' \item Mark data for exclusion according to the following criteria:
#' \itemize{
#'     \item Exclude any remaining observations where both Value and LOQ are NA
#'     \item For any cases where `N_Subjects` is NA, impute `N_Subjects` = 1
#'     \item Exclude any remaining observations with `N_Subjects` > 1 and `Value_SD` still NA. (This should never occur, if SD imputation is performed, but just in case.)
#'     \item Exclude any observations with `N_Subjects` > 1 where reported `Value` is NA, because log-likelihood for non-detect multi-subject observations has not been implemented.
#'     \item Exclude any observations with NA `Time` values
#'     \item Exclude any observations with `Dose` = 0
#'     }
#' \item Apply any time transformations specified by user
#' \item Scale concentration by `ratio_conc_dose`
#' \item Apply any concentration transformations specified by the user.
#' \item If `Series_ID` is not included, then assign it as NA
#' \item Create variable `pLOQ` and set it equal to `LOQ`
#' }
#'
#' @param obj A `pk` object
#' @param ... Additional arguments. Not in use currently.
#' @return The same `pk` object, with added elements `data` (containing the
#'   cleaned, gap-filled data) and `data_info` (containing summary information
#'   about the data, e.g. number of observations by route, media,
#'   detect/nondetect; empirical tmax, time of peak concentration for oral data;
#'   number of observations before and after empirical tmax)
#' @author John Wambaugh, Caroline Ring, Christopher Cook, Gilberto Padilla Mercado
#' @export
do_preprocess.pk <- function(obj, ...) {
  suppress_messages <- obj$settings_preprocess$suppress.messages
  if (!suppress_messages) {
    cli_inform("do_preprocess.pk(): Pre-processing data")
  }

  objname <- deparse(substitute(obj))
  if (obj$status >= status_preprocess) {
    cli_warn(c(
        "{objname} current status is {obj$status}.",
        "i" = paste(
          "{.fn do_preprocess.pk} will reset its status to {status_preprocess}",
          " and any results from later workflow stages will be lost."
        )
    ))
  }

  if (is.null(obj$data_original)) {
    cli_warn(c("do_preprocess.pk(): Original data is NULL",
               "i" = "Returning {objname}."))
    obj$data <- NULL
    obj$data_info <- NULL
    obj$status <- status_preprocess
    return(obj)
  } else {
    # coerce to data.frame (in case it is a tibble or data.table or whatever)
    data_original <- as.data.frame(obj$data_original)

    # rename variables using the mapping
    data <- as.data.frame(
      sapply(obj$mapping,
        rlang::eval_tidy,
        data_original,
        simplify = FALSE,
        USE.NAMES = TRUE)
    )

    ### Coerce Species, Route, and Media to lowercasei
    if (!suppress_messages) {
      cli_inform(
        "Species, Route, and Media will be coerced to lowercase."
      )
    }
    # Check hierarchical grouping structure
    check_group_hierarchy(obj)

    data$Species <- tolower(data$Species)
    data$Route <- tolower(data$Route)
    data$Media <- tolower(data$Media)

    # Check to make sure the data include only route_keep and media_keep
    routes <- unique(data$Route)
    media <- unique(data$Media)

    if (!(
      all(routes %in% obj$settings_preprocess$routes_keep) &&
      all(media %in% obj$settings_preprocess$media_keep)
    )) {
      cli_abort(c(
        "do_preprocess.pk(): data contains unsupported media and/or routes.",
        "i" = "Supported media: {obj$settings_preprocess$media_keep}",
        "Media in data: {media}",
        "i" = "Supported route: {obj$settings_preprocess$routes_keep}",
        "Route{?s} in data: {route}"

      ))
    }

    # Check to make sure there is only one value for each set of units
    time_units_ <- unique(data$Time.Units) # Underscore to prevent name conflict with internal data
    value_units <- unique(data$Value.Units)
    weight_units <- unique(data$Weight.Units)
    dose_units <- unique(data$Dose.Units)

    # Externalized this for semantic error messages
    non_unique_units <- unname(
      lengths(list(time_units_, value_units, weight_units, dose_units)) > 1
      )

    if (any(non_unique_units)) {
      non_unique_units <- c("Time.Units",
                            "Value.Units",
                            "Weight.Units",
                            "Dose.Units")[non_unique_units]
      cli_abort(c(
        "do_preprocess.pk(): data contains multiple units for {non_unique_units}",
        "Please ensure that {non_unique_units} {?has/have} unique units."
      ))
    }

    # Check to make sure that time units are in allowable list
    if (!all(unique(time_units_ %in% time_units))) {
      cli_abort(c(
        "do_preprocess.pk(): Data has Time.Units {time_units_}",
        "i" = "Please ensure that Time.Units is one of the allowable units: {time_units}",
        " " = "(see built-in data object {.envvar time_units})."
      ))
    }

    # If data has passed all these initial checks, then proceed with pre-processing
    data_group <- obj$data_group
    summary_group <- obj$settings_data_info$summary_group

    data_group_chr <- sapply(data_group, rlang::as_label)

    n_grps <- dplyr::group_by(data, !!!data_group) |>
      dplyr::group_keys() |>
      dplyr::n_distinct()

    if (!suppress_messages) {
      ### display messages describing loaded data
      cli_inform(c(
        "{nrow(data)} concentration vs time observation{?s} loaded.",
        paste(
          "Number of unique data groups",
          "(unique combinations of {data_group_chr}): {n_grps}"
        )
      ))
    }

    ### Coerce all 'Value' values to be numeric and say so
    if (!is.numeric(data$Value)) {
      value_num <- as.numeric(data$Value)
      old_na <- sum(is.na(data$Value) | !nzchar(data$Value))
      new_na <- sum(is.na(value_num))
      if (!suppress_messages) {
        cli::cli_inform(c(
          paste(
            "Column {.envvar Value} converted from",
            "{.cls {class(data$Value)}}",
            "to {.cls {class(value_num)}}."
          ),
          "i" = "Pre-conversion NAs and blanks: {old_na}",
          "i" = "Post-conversion NAs and blanks: {new_na}"
        ))
      }
      data$Value <- value_num
      rm(value_num, old_na, new_na)
    }

    ### Coerce all 'Value_SD' values to be numeric and say so
    if (!is.numeric(data$Value_SD)) {
      valuesd_num <- as.numeric(data$Value_SD)
      old_na <- sum(is.na(data$Value_SD) | !nzchar(data$Value_SD))
      new_na <- sum(is.na(valuesd_num))
      if (!suppress_messages) {
        cli::cli_inform(c(
          paste(
            "Column {.envvar Value_SD} converted from",
            "{.cls {class(data$Value_SD)}}",
            "to {.cls {class(valuesd_num)}}."
          ),
          "i" = "Pre-conversion NAs and blanks: {old_na}",
          "i" = "Post-conversion NAs and blanks: {new_na}"
        ))
      }
      data$Value_SD <- valuesd_num
      rm(valuesd_num, old_na, new_na)
    }

    ### Coerce all 'LOQ' values to be numeric and say so
    if (!is.numeric(data$LOQ)) {
      loq_num <- as.numeric(data$LOQ)
      old_na <- sum(is.na(data$LOQ) | !nzchar(data$LOQ))
      new_na <- sum(is.na(loq_num))
      if (!suppress_messages) {
        cli::cli_inform(c(
          paste(
            "Column {.envvar LOQ} converted from",
            "{.cls {class(data$LOQ)}}",
            "to {.cls {class(loq_num)}}."
          ),
          "i" = "Pre-conversion NAs and blanks: {old_na}",
          "i" = "Post-conversion NAs and blanks: {new_na}"
        ))
      }
      data$LOQ <- loq_num
      rm(loq_num, old_na, new_na)
    }


    ### coerce 'Dose' values to numeric and say so
    if (!is.numeric(data$Dose)) {
      dose_num <- as.numeric(data$Dose)
      old_na <- sum(is.na(data$Dose) | !nzchar(data$Dose))
      new_na <- sum(is.na(dose_num))
      if (!suppress_messages) {
        cli::cli_inform(c(
          paste(
            "Column {.envvar Dose} converted from",
            "{.cls {class(data$Dose)}}",
            "to {.cls {class(dose_num)}}."
          ),
          "i" = "Pre-conversion NAs and blanks: {old_na}",
          "i" = "Post-conversion NAs and blanks: {new_na}"
        ))
      }
      data$Dose <- dose_num
      rm(dose_num, old_na, new_na)
    }

    ### coerce 'Time' values to numeric and say so
    if (!is.numeric(data$Time)) {
      time_num <- as.numeric(data$Time)
      old_na <- sum(is.na(data$Time) | !nzchar(data$Time))
      new_na <- sum(is.na(time_num))
      if (!suppress_messages) {
        cli::cli_inform(c(
          paste(
            "Column {.envvar Time} converted from",
            "{.cls {class(data$Time)}}",
            "to {.cls {class(time_num)}}."
          ),
          "i" = "Pre-conversion NAs and blanks: {old_na}",
          "i" = "Post-conversion NAs and blanks: {new_na}"
        ))
      }
      data$Time <- time_num
      rm(time_num, old_na, new_na)
    }

    ### coerce 'N_Subjects' values to numeric and say so
    if (!is.numeric(data$N_Subjects)) {
      N_Subjects_num <- as.numeric(data$N_Subjects)
      old_na <- sum(is.na(data$N_Subjects) |
                      !nzchar(data$N_Subjects))
      new_na <- sum(is.na(N_Subjects_num))
      if (!suppress_messages) {
        cli::cli_inform(c(
          paste(
            "Column {.envvar N_Subjects} converted from",
            "{.cls {class(data$N_Subjects)}}",
            "to {.cls {class(N_Subjects_num)}}."
          ),
          "i" = "Pre-conversion NAs and blanks: {old_na}",
          "i" = "Post-conversion NAs and blanks: {new_na}"
        ))
      }
      data$N_Subjects <- N_Subjects_num
      rm(N_Subjects_num, old_na, new_na)
    }

    data$Value_orig <- data$Value
    data$LOQ_orig <- data$LOQ
    data$Value_SD_orig <- data$Value_SD

    # Coerce any negative Value to NA
    data[(data$Value < 0) %in% TRUE, "Value"] <- NA_real_
    if (!suppress_messages && anyNA(data$Value)) {
        cli_inform(
          c(
            "If Value < 0, replacing Value with NA. ",
            "i" = "{sum(is.na(data$Value))} Values were be replaced with NA.\n"
          )
        )
    }

    # Coerce any negative LOQ to NA
    data[(data$LOQ < 0) %in% TRUE, "LOQ"] <- NA_real_
    if (!suppress_messages && anyNA(data$LOQ)) {
        cli_inform(c(
            "If LOQ < 0, replacing LOQ with NA.",
            "i" = "{sum(is.na(data$LOQ))} LOQs were be replaced with NA."
        ))
    }

    # Coerce any non-positive Value_SD to NA
    data[(data$Value_SD < 0) %in% TRUE, "Value_SD"] <- NA_real_
    if (!suppress_messages && anyNA(data$Value_SD)) {
        cli_inform(c(
            "If Value_SD <= 0, replacing Value_SD with NA.",
            "i" = "{sum(is.na(data$Value_SD))} Value_SDs were be replaced with NA."
        ))
    }

    # Coerce any negative Time to NA
    data[(data$Time < 0) %in% TRUE, "Time"] <- NA_real_
    if (!suppress_messages && anyNA(data$Time)) {
        cli_inform(c(
            "If Time < 0, replacing Time with NA.",
            "i" = "{sum(is.na(data$Time))} Times were be replaced with NA."
        ))
    }

    # If any non-NA Value is currently less than its non-NA LOQ, replace it with NA
    if (!suppress_messages && any((data$Value <= data$LOQ) %in% TRUE)) {
        cli_inform(c(
            "If Value <= LOQ, replacing Value with NA.",
            "i" = "{sum((data$Value <= data$LOQ) %in% TRUE)} Values were be replaced with NA."
        ))
    }
    data$Value <- ifelse((data$Value <= data$LOQ) %in% TRUE, NA_real_, data$Value)

    # Impute LOQ:
    # as calc_loq_factor * minimum non-NA value in each loq_group
    if (obj$settings_preprocess$impute_loq %in% TRUE) {
      if (anyNA(data$LOQ) && !suppress_messages) {
        cli_inform(c(
          paste("Estimating {sum(is.na(data$LOQ))} missing LOQ{?s} as ",
                "{obj$settings_preprocess$calc_loq_factor}",
                "* minimum non-NA, non-zero Value for each unique combination of ",
                "{sapply(obj$settings_preprocess$loq_group, as_label)}"
          )
        ))
      }
      data <- dplyr::group_by(data, !!!obj$settings_preprocess$loq_group) |>
        dplyr::mutate(LOQ_orig = LOQ,
                      LOQ = ifelse(is.na(LOQ_orig), {
                        if (any((Value > 0) %in% TRUE)) {
                          min(Value[(Value > 0) %in% TRUE], na.rm = TRUE) *
                            obj$settings_preprocess$calc_loq_factor
                        } else {
                          NA_real_
                        }
                      }, LOQ_orig)) |>
        dplyr::ungroup() |>
        as.data.frame()
    } # end if impute_loq %in% TRUE

    if (!suppress_messages && any((data$Value < data$LOQ) %in% TRUE)) {
        cli_inform(c(
            "Converting 'Value' values of less than LOQ to NA.",
            "i" = "{sum((data$Value < data$LOQ) %in% TRUE)} value{?s} will be converted."
        ))
    }
    data[(data$Value < data$LOQ) %in% TRUE, "Value"] <- NA_real_

    data$exclude <- FALSE
    data$exclude_reason <- NA_character_

    # Exclude any remaining cases where both Value and LOQ are NA
    # because if we don't have a value or an LOQ, we can't do anything
    if (anyNA(data$Value) && anyNA(data$LOQ)) {

      if (!suppress_messages) {
        cli_inform(c(
            "Excluding observations where both Value and LOQ were NA.",
            "{sum(is.na(data$Value) & is.na(data$LOQ))}, observation{?s} were excluded.\n"
        ))
      }
      data$exclude <- ifelse(is.na(data$Value) & is.na(data$LOQ), TRUE, data$exclude)

      data$exclude_reason <- ifelse(
        is.na(data$Value) & is.na(data$LOQ),
        paste2(data$exclude_reason, "Value & LOQ both NA", sep = "; "),
        data$exclude_reason
      )
    }

    # For any cases where N_Subjects is NA, impute N_Subjects = 1
    if (anyNA(data$N_Subjects)) {
      data$N_Subjects_orig <- data$N_Subjects
      if (!suppress_messages) {
        cli_inform(c(
          "N_Subjects is NA for {sum(is.na(data$N_Subjects))} observation{?s}.",
          "i" = "These will be assumed = 1."
        ))
      }
      data[is.na(data$N_Subjects), "N_Subjects"] <- 1
    }

    # for anything with N_Subjects == 1, set Value_SD to 0
    if (any(data$N_Subjects == 1)) {
      if (!suppress_messages) {
        cli_inform(c(
            "N_Subjects is 1 for {sum(data$N_Subjects == 1)} observation{?s}.",
            "i" = "Value_SD will be set to 0 for these observations.\n"
        ))
      }
      data[data$N_Subjects == 1, "Value_SD"] <- 0
    }

    # Impute missing SDs
    data$Value_SD_orig <- data$Value_SD
    if (obj$settings_preprocess$impute_sd %in% TRUE) {
      if (any((data$N_Subjects > 1) %in% TRUE & is.na(data$Value_SD)) && !suppress_messages) {
        # number of SDs to be estimated
        n_sd_est <- sum((data$N_Subjects > 1) %in% TRUE &
                          is.na(data$Value_SD))
        cli_inform(c(
          paste(
            "Estimating missing concentration SDs ",
            "for data points with N_Subjects > 1 as ",
            "minimum non-missing SD for each group of data",
            "given by the unique combination of variables",
            "{sapply(obj$settings_preprocess$sd_group, as_label)}",
            ". If all SDs are missing in a group, ",
            "SD for each observation will be imputed as",
            " 0. "
          ),
          "i" = "{n_sd_est} missing SD{?s} will be estimated."
        ))
      }
      data <- dplyr::group_by(data, !!!obj$settings_preprocess$sd_group) |>
        dplyr::mutate(
          Value_SD_orig = Value_SD,
          Value_SD = ifelse(
            is.na(Value_SD_orig) &
              N_Subjects > 1,
            ifelse(
              rep(all(is.na(
                Value_SD_orig
              )), dplyr::n()),
              0,
              min(Value_SD_orig, na.rm = TRUE)
            ),
            Value_SD_orig
          )
        ) |>
        dplyr::ungroup() |>
        as.data.frame()
    } # end if impute_sd %in% TRUE

# Exclude any remaining multi-subject observations where SD is NA
if (any((data$N_Subjects > 1) %in% TRUE & is.na(data$Value_SD))) {
  n_exclude <- sum((data$N_Subjects > 1) %in% TRUE & is.na(data$Value_SD))

  data$exclude <- ifelse((data$N_Subjects > 1) %in% TRUE &
                           is.na(data$Value_SD),
                         TRUE,
                             data$exclude)
      data$exclude_reason <- ifelse((data$N_Subjects > 1) %in% TRUE &
                                      is.na(data$Value_SD),
                                    paste2(
                                      data$exclude_reason,
                                      "N_Subjects > 1 and Value_SD is NA",
                                      sep = "; "
                                    ),
                                    data$exclude_reason)
      if (!suppress_messages) {
        cli_inform(c(
          paste(
            "Excluding remaining observations with N_Subjects > 1 ",
            "where reported SD is NA. "
          ),
          "i" = "{sum(exclude)} observation{?s} were excluded."
        ))
      }
    }

    if (!suppress_messages) {
      n_grps <- dplyr::group_by(data, !!!data_group) |>
        dplyr::group_keys() |>
        dplyr::n_distinct()
      ### display messages describing loaded data
      cli_inform(c(
        "Remaining observations: {nrow(data)}",
        paste(
          "Number of unique data groups",
          "(unique combinations of {data_group_chr}): {n_grps}"
        )
      ))
    }
    # Update the remaining observations
    obs_num <- nrow(data)

    # Exclude any remaining multi-subject observations where Value is NA
    if (any((data$N_Subjects > 1) %in% TRUE & is.na(data$Value))) {
      n_exclude <- sum((data$N_Subjects > 1) %in% TRUE & is.na(data$Value))

      data$exclude <- ifelse((data$N_Subjects > 1) %in% TRUE &
                               is.na(data$Value),
                             TRUE,
                             data$exclude)

      data$exclude_reason <- ifelse((data$N_Subjects > 1) %in% TRUE &
                                      is.na(data$Value),
                                    paste2(data$exclude_reason,
                                           "N_Subjects > 1 and Value is NA",
                                           sep = "; "),
                                    data$exclude_reason
      )

      if (!suppress_messages) {
        cli_inform(c(paste(
          "Excluding observations with N_Subjects > 1 ",
          "where reported Value is NA ",
          "(because log-likelihood for non-detect multi-subject observations ",
          "has not been implemented). "
        ),
          "i" = "{n_exclude} observation{?s} will be excluded."
        ))
      }

    }

    if (!suppress_messages && (nrow(data) != obs_num)) {
      n_grps <- dplyr::group_by(data, !!!data_group) |>
        dplyr::group_keys() |>
        dplyr::n_distinct()
      ### display messages describing loaded data
      cli_inform(c(
        "Remaining observations: {nrow(data)}",
        paste(
          "Number of unique data groups",
          "(unique combinations of {data_group_chr}): {n_grps}"
        )
      ))
      obs_num <- nrow(data)
    }

    # Exclude any NA time values
    if (anyNA(data$Time)) {
      data$exclude <- ifelse(is.na(data$Time), TRUE, data$exclude)

      data$exclude_reason <- ifelse(
        is.na(data$Time),
        paste2(data$exclude_reason, "Time is NA", sep = "; "),
        data$exclude_reason
      )
      if (!suppress_messages) {
        cli_inform(c(
            "Excluding observations with NA time values.",
            "i" = "{sum(is.na(data$Time))} observation{?s} will be excluded."
        ))
      }

    }

    if (!suppress_messages && (nrow(data) != obs_num)) {
      n_grps <- dplyr::group_by(data, !!!data_group) |>
        dplyr::group_keys() |>
        dplyr::n_distinct()
      ### display messages describing loaded data
      cli_inform(c(
        "Remaining observations: {nrow(data)}",
        paste(
          "Number of unique data groups",
          "(unique combinations of {data_group_chr}): {n_grps}"
        )
      ))
      obs_num <- nrow(data)
    }

    # Exclude any Dose = 0 observations
    if (any(data$Dose <= .Machine$double.eps)) {

      data$exclude <- ifelse(data$Dose < .Machine$double.eps, TRUE, data$exclude)

      data$exclude_reason <- ifelse(
        data$Dose < .Machine$double.eps,
        paste2(data$exclude_reason, "Dose == 0", sep = "; "),
        data$exclude_reason
      )
      if (!suppress_messages) {
        cli_inform(c(
          "Excluding observations with Dose == 0 (control observations). ",
          "i" = "{sum(data$Dose <= .Machine$double.eps)} observation{?s} will be excluded."
        ))
      }

    }

    if (!suppress_messages && (nrow(data) != obs_num)) {
      n_grps <- dplyr::group_by(data, !!!data_group) |>
        dplyr::group_keys() |>
        dplyr::n_distinct()
      ### display messages describing loaded data
      cli_inform(c(
        "Remaining observations: {nrow(data)}",
        paste(
          "Number of unique data groups",
          "(unique combinations of {data_group_chr}): {n_grps}"
        )
      ))
      obs_num <- nrow(data)
    }

    # apply time transformation
    data$Time_orig <- data$Time
    data$Time.Units_orig <- data$Time.Units

    # first, default to identity transformation if none is specified
    if (is.null(obj$scales$time$new_units)) {
      obj$scales$time$new_units <- "identity"
    }


    from_units <- unique(data$Time.Units)
    to_units <- ifelse(
      obj$scales$time$new_units %in% "identity",
      from_units,
      obj$scales$time$new_units
    )


    # Needs to be grouped
    if (obj$scales$time$new_units %in% "auto") {
      to_units <- data |>
        dplyr::group_by(!!!data_group) |>
        dplyr::mutate(AUTO_UNITS = auto_units(y = Time, from = Time.Units_orig)) |>
        dplyr::pull(AUTO_UNITS)
    }

    # Needs to be iterated
    if (all(unique(to_units) %in% from_units)) {
      data$Time_trans <- data$Time

    } else {
      if (!suppress_messages) {
        for (convert_message in as.numeric(which(!(unique(to_units) %in% from_units)))) {
          cli_inform(
            "Converting time from {from_units} to {unique(to_units)[convert_message]}"
          )
        }

      }
      data$NEW_UNITS <- to_units
      data <- data |>
        dplyr::rowwise() |>
        dplyr::mutate(Time_trans = tryCatch(
          convert_time(
            x = Time,
            from = from_units,
            to = NEW_UNITS,
            inverse = FALSE
          ),
          error = function(err) {
            warning(
                "invivopkfit::do_preprocess.pk():",
                "Error in transforming time using convert_time():",
                err$message
            )
            return(NA_real_)
          }
        )) |>
        dplyr::select(-NEW_UNITS)
    }
    data$Time_trans.Units <- to_units
    # End grouping


    # Scale concentration by ratio_conc_dose

    if (is.null(obj$scales$conc$expr)) {
      obj$settings_preprocess$ratio_conc_dose <- 1
    }

    ratio_conc_dose <- obj$settings_preprocess$ratio_conc_dose

    if (!suppress_messages && !(ratio_conc_dose %in% 1)) {
      cli_inform(c(paste(
        "Concentrations, LOQs, and concentration SDs are being scaled by ",
        "ratio_conc_dose = {ratio_conc_dose}."
      )))
    }


    if ("Conc" %in% names(data) && !suppress_messages) {
      cli_inform(c(
        "Harmonized variable `Conc` already exists in data!. ",
        "i" = paste(
          "It will be overwritten as ",
          "`pmax(Value/ratio_conc_dose, LOQ/ratio_conc_dose, na.rm = TRUE)`."
        )
      ))
    }
    data$Conc <- pmax(data$Value / ratio_conc_dose,
                      data$LOQ / ratio_conc_dose,
                      na.rm = TRUE)

    if ("Detect" %in% names(data)) {
      if (!suppress_messages) {
        cli_inform(c(
            "Harmonized variable `Detect`already exists in data! ",
            "i" = "It will be overwritten as `is.na(Value)`."
        ))
      }
    }
    data$Detect <- !is.na(data$Value)

    if ("Conc_SD" %in% names(data) && !suppress_messages) {
        cli_inform(c(
            "Harmonized variable `Conc_SD` already exists in data! ",
            "i" = paste(
              "It will be overwritten as {.code Value_SD/ratio_conc_dose}",
              "if {.envvar Detect} is {.code TRUE}",
              "otherwise it is set to {.code NA}."
            )
        ))
    }
    data$Conc_SD <- data$Value_SD / ratio_conc_dose

    if ("Conc.Units" %in% names(data) && !suppress_messages) {
        cli_inform(c(
            "Harmonized variable `Conc.Units` already exists in data!",
            "i" = "It will be overwritten as {.code Value.Units/ratio_conc_dose}."
        ))
    }
    data$Conc.Units <- ifelse(
      ratio_conc_dose == 1,
      data$Value.Units,
      paste0(data$Value.Units, "/", ratio_conc_dose)
    )

    # where it is excreta make cumulative sum, if needed
    if ("excreta" %in% data$Media && !suppress_messages) {
      cli_inform(c(
        "{.emph excreta} is a value {.envvar Media}",
        "i" = paste(
          "Changing all {.envvar Conc} values per subject to cummulative concentrations",
          "where {.envvar Media} is excreta"
        )
      ))
    }
    data <- data |>
      dplyr::group_by(!!!summary_group, Subject_ID) |>
      dplyr::arrange(Time) |>
      dplyr::mutate(
        Conc = ifelse(
          Media %in% "excreta",
          cumsum(Conc),
          Conc
        )
      )


    # apply concentration transformation.
    #
    # .conc is a placeholder that refers to any concentration variable (Conc,
    # value, Value_SD, Conc_SD, LOQ).

    # apply conc transformation
    # use tidy evaluation: Specify an expression to evaluate in the context of a data.frame
    # Transform Conc (this is either the measured value or the LOQ, depending on Detect)
    # If no conc transformation specified, assume identity
    if (is.null(obj$scales$conc$expr)) {
      obj$scales$conc$dose_norm <- FALSE
      obj$scales$conc$log10_trans <- FALSE
      obj$scales$conc$expr <- rlang::new_quosure(quote(.conc), env = caller_env())
    }

    if (!suppress_messages) {
      cli_inform(c(
          "Applying transformations to concentration variables:",
          "dose_norm = {obj$scales$conc$dose_norm}",
          "log10_trans = {obj$scales$conc$log10_trans}"
      ))
      if ("Conc_trans" %in% names(data)) {
        cli_inform(c(
          "{.emph Warning}: variable `Conc_trans` already exists in the data. ",
          "It will be overwritten!"
        ))
      }
      if ("Conc_SD_trans" %in% names(data)) {
        cli_inform(c(
          "{.emph Warning}: variable `Conc_SD_trans` already exists in the data. ",
          "It will be overwritten!"
        ))
      }
      if ("Conc_trans.Units" %in% names(data)) {
        cli_inform(c(
          "{.emph Warning}: variable `Conc_trans.Units` already exists in the data.",
          "It will be overwritten!"
        ))
      }
    }

    data$Conc_trans <- rlang::eval_tidy(
      expr = obj$scales$conc$expr,
      data = cbind(data, data.frame(.conc = data$Conc)))

    # Transform Conc_SD
    data$Conc_SD_trans <- rlang::eval_tidy(
      obj$scales$conc$expr,
      data = cbind(data, data.frame(.conc = data$Conc_SD)))

    # Record new conc units
    data$Conc_trans.Units <- gsub(
      x = gsub(
        x = rlang::as_label(obj$scales$conc$expr),
        pattern = ".conc",
        replacement = paste0("(", unique(data$Conc.Units), ")"),
        fixed = TRUE
      ),
      pattern = "Dose",
      replacement = paste0("(", unique(data$Dose.Units), ")"),
      fixed = TRUE
    )

    # If Series ID is not included, then assign it as NA
    if (!("Series_ID" %in% names(data))) {
      data$Series_ID <- NA_integer_
    }

    # Make sure that Reference is a list of characters
    data$Reference <- as.character(data$Reference)

    # Create prediction LOQ column
    data$pLOQ <- data$LOQ


    # add data & data info to object AND a group id column
    data <- dplyr::ungroup(data)
    obj$data <- tidyr::unite(data, "GRP_ID", !!!data_group, remove = FALSE)

    if (is.null(obj$settings_preprocess$keep_data_original)) {
      obj$settings_preprocess$keep_data_original <- TRUE
    }

    if (obj$settings_preprocess$keep_data_original == FALSE) {
      obj$data_original <- NULL
    }

    obj$status <- status_preprocess # preprocessing complete

    return(obj)
  }
}

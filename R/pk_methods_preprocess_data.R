#' Pre-process data
#'
#' Pre-process data for a `pk` object
#'
#' Data pre-processing for an object `obj` includes the following steps, in order:
#' 1. Coerce data to class `data.frame` (if it is not already)
#' 2. Rename variables to harmonized "`invivopkfit` aesthetic" variable names, using `obj$mapping`
#' 3. Check that the data includes only one chemical and one species.
#' 4. Check that the data includes only routes in `obj$settings_preprocess$routes_keep` and media in `obj$settings_preprocess$media_keep`
#' 5. Check that the data includes only one unit for concentration, one unit for time, and one unit for dose.
#' 6. Coerce Value, Value_SD, LOQ, Dose, and Time to numeric, if they are not already.
#' 7. Coerce Species, Route, and Media to lowercase.
#' 8. Replace any negative Value, Value_SD, Dose, or Time with `NA`
#' 9. Replace any negative
#'
#'
#' @param obj A `pk` object
#' @return The same `pk` object, with added elements `data` (containing the
#'   cleaned, gap-filled data) and `data_info` (containing summary information
#'   about the data, e.g. number of observations by route, media,
#'   detect/nondetect; empirical tmax, time of peak concentration for oral data;
#'   number of observations before and after empirical tmax)
#' @author Caroline Ring
#' @importFrom magrittr `%>%`
#' @export
preprocess_data.pk <- function(obj){

  if(!obj$settings_preprocess$suppress.messages){
    message("preprocess_data.pk(): Pre-processing data\n")
  }

  objname <- deparse(substitute(obj))
  if(obj$status >= status_preprocess){
    warning(paste0(objname,
                  " current status is ",
                  obj$status,
                  ". preprocess_data() will reset its status to ",
                  status_preprocess,
                  ". Any results from later workflow stages will be lost.\n"))
  }

  if(is.null(obj$data_original)){
    message("preprocess_data.pk(): Original data is NULL")
    obj$data <- NULL
    obj$data_info <- NULL
    obj$status <- status_preprocess
    return(obj)
  }else{

    #coerce to data.frame (in case it is a tibble or data.table or whatever)
    data_original <- as.data.frame(obj$data_original)

    #rename variables using the mapping
    data <- as.data.frame(sapply(obj$mapping,
                                 function(x) rlang::eval_tidy(x, data_original),
                                 simplify = FALSE,
                                 USE.NAMES = TRUE)
    )

    ### Coerce Species, Route, and Media to lowercase
    if(!obj$settings_preprocess$suppress.messages){
      message(paste("Species, Route, and Media",
                    "will be coerced to lowercase.\n"))
    }
    data$Species <- tolower(data$Species)
    data$Route <- tolower(data$Route)
    data$Media <- tolower(data$Media)

    # #Check to make sure the data include only one Chemical and Species. Stop
    # #with an error otherwise.
    chems <- unique(data$Chemical)
    species <- unique(data$Species)
    #
    # nchem <- length(chems)
    # nspecies <- length(species)
    #
    # if(!(nchem %in% 1 & nspecies %in% 1)){
    #   stop(paste("preprocess_data.pk(): data contains multiple chemicals and/or multiple species.",
    #              "Unique chemicals in this data:",
    #              paste(chems, collapse = "; "),
    #              "Unique Species in this data:",
    #              paste(species, collapse = "; "),
    #              sep = "\n"))
    # }

    #Check to make sure the data include only route_keep and media_keep
    routes <- unique(data$Route)
    media <- unique(data$Media)

    if( !(all(routes %in% obj$settings_preprocess$routes_keep) &
          all(media %in% obj$settings_preprocess$media_keep)) ){
      stop(paste("preprocess_data.pk(): data contains unsupported media and/or routes.",
                 paste("Supported media:", paste(obj$settings_preprocess$media_keep, collapse = "; ")),
                 paste("Media in data:", paste(media, collapse = "; ")),
                 paste("Supported routes:", paste(obj$settings_preprocess$routes_keep, collapse = "; ")),
                 paste("Routes in data:", paste(routes, collapse = "; ")),
                 sep = "\n"))
    }




    #Check to make sure there is only one value for each set of units
    time_units <- unique(data$Time.Units)
    value_units <- unique(data$Value.Units)
    weight_units <- unique(data$Weight.Units)
    dose_units <- unique(data$Dose.Units)

    if(any(sapply(list(time_units,
                       value_units,
                       weight_units,
                       dose_units),
                  function(x) length(x) > 1)
    )){
      stop(paste("preprocess_data.pk(): data contains multiple units for one or more variables.",
                 paste("Time units:", paste(time_units, collapse = "; "), sep = " ",),
                 paste("Concentration units:", paste(value_units, collapse = "; "), sep = " "),
                 paste("Weight units:", paste(weight_units, collapse = "; "), sep = " "),
                 paste("Dose units:", paste(dose_units, collapse = "; "), sep = " "),
                 sep = "\n"))
    }

    #Check to make sure that time units are in allowable list
    if(!all(unique(data$Time.Units %in% time_units))){
      stop(paste0("preprocess_data.pk(): Data has Time.Units ",
                 unique(data$Time.Units),
                 " that are not on the list of allowable time units (see built-in data object `time_units`):\n",
                 paste(time_units, collapse = "\n")))
    }

    # If data has passed all these initial checks, then proceed with pre-processing

    if(!obj$settings_preprocess$suppress.messages){
      ### display messages describing loaded data
      message(
        paste0(
          paste(
          paste(nrow(data),
                "concentration vs. time observations loaded."),
          paste("Chemicals:", chems),
          paste("Species:", species),
          paste("Routes:", routes),
          paste("Media:", media),
          sep = "\n"
        ),
        "\n")
      )
    }

    ### Coerce all 'Value' values to be numeric and say so
    if (!is.numeric(data$Value))
    {
      value_num <- as.numeric(data$Value)
      old_na <- sum(is.na(data$Value) | !nzchar(data$Value))
      new_na <- sum(is.na(value_num))
      if(!obj$settings_preprocess$suppress.messages){
        message(paste0("Column \"Value\" converted from ",
                       class(data$Value),
                       " to numeric. ",
                       "Pre-conversion NAs and blanks: ",
                       old_na,
                       ". Post-conversion NAs: ",
                       new_na, ".",
                       "\n"))
      }
      data$Value <- value_num
      rm(value_num, old_na, new_na)
    }

    ### Coerce all 'Value_SD' values to be numeric and say so
    if (!is.numeric(data$Value_SD))
    {
      valuesd_num <- as.numeric(data$Value_SD)
      old_na <- sum(is.na(data$Value_SD) | !nzchar(data$Value_SD))
      new_na <- sum(is.na(valuesd_num))
      if(!obj$settings_preprocess$suppress.messages){
        message(paste0("Column \"Value_SD\" converted from ",
                       class(data$Value),
                       " to numeric. ",
                       "Pre-conversion NAs and blanks: ",
                       old_na,
                       ". Post-conversion NAs: ",
                       new_na, ".\n"))
      }
      data$Value_SD <- valuesd_num
      rm(valuesd_num, old_na, new_na)
    }

    ### Coerce all 'LOQ' values to be numeric and say so
    if (!is.numeric(data$LOQ))
    {
      loq_num <- as.numeric(data$LOQ)
      old_na <- sum(is.na(data$LOQ) | !nzchar(data$LOQ))
      new_na <- sum(is.na(loq_num))
      if(!obj$settings_preprocess$suppress.messages){
        message(paste0("Column \"LOQ\" converted from ",
                       class(data$LOQ),
                       " to numeric. ",
                       "Pre-conversion NAs and blanks: ",
                       old_na,
                       ". Post-conversion NAs: ",
                       new_na, ".\n"))
      }
      data$LOQ <- loq_num
      rm(loq_num, old_na, new_na)
    }


    ### coerce 'Dose' values to numeric and say so
    if (!is.numeric(data$Dose))
    {
      dose_num <- as.numeric(data$Dose)
      old_na <- sum(is.na(data$Dose) | !nzchar(data$Dose))
      new_na <- sum(is.na(dose_num))
      if(!obj$settings_preprocess$suppress.messages){
        message(paste0("Column \"Dose\" converted from ",
                       class(data$Dose),
                       " to numeric. ",
                       "Pre-conversion NAs and blanks: ",
                       old_na,
                       ". Post-conversion NAs: ",
                       new_na, ".\n"))
      }
      data$Dose <- dose_num
      rm(dose_num, old_na, new_na)
    }

    ### coerce 'Time' values to numeric and say so
    if (!is.numeric(data$Time))
    {
      time_num <- as.numeric(data$Time)
      old_na <- sum(is.na(data$Time) | !nzchar(data$Time))
      new_na <- sum(is.na(time_num))
      if(!obj$settings_preprocess$suppress.messages){
        message(paste0("Column \"Time\" converted from ",
                       class(data$Time),
                       " to numeric. ",
                       "Pre-conversion NAs and blanks: ",
                       old_na,
                       ". Post-conversion NAs: ",
                       new_na, ".\n"))
      }
      data$Time <- time_num
      rm(time_num, old_na, new_na)
    }

    ### coerce 'N_Subjects' values to numeric and say so
    if (!is.numeric(data$N_Subjects))
    {
      N_Subjects_num <- as.numeric(data$N_Subjects)
      old_na <- sum(is.na(data$N_Subjects) | !nzchar(data$N_Subjects))
      new_na <- sum(is.na(N_Subjects_num))
      if(!obj$settings_preprocess$suppress.messages){
        message(paste0("Column \"N_Subjects\" converted from ",
                       class(data$N_Subjects),
                       " to numeric. ",
                       "Pre-conversion NAs and blanks: ",
                       old_na,
                       ". Post-conversion NAs: ",
                       new_na, ".\n"))
      }
      data$N_Subjects <- N_Subjects_num
      rm(N_Subjects_num, old_na, new_na)
    }



data$Value_orig <- data$Value
data$LOQ_orig <- data$LOQ
data$Value_SD_orig <- data$Value_SD

    #Coerce any negative Value to NA
    if(!obj$settings_preprocess$suppress.messages){
      if(any((data$Value<0) %in% TRUE)){
        message(paste('If value < 0, replacing Value with NA.',
                      sum((data$Value < 0) %in% TRUE),
                      "Values will be replaced with NA.\n"))
      }
    }
    data[(data$Value<0) %in% TRUE, "Value"] <- NA_real_

    #Coerce any negative LOQ to NA
    if(!obj$settings_preprocess$suppress.messages){
      if(any((data$LOQ<0) %in% TRUE)){
        message(paste('If LOQ < 0, replacing LOQ with NA.',
                      sum((data$LOQ < 0) %in% TRUE),
                      "LOQs will be replaced with NA.\n"))
      }
    }
    data[(data$LOQ<0) %in% TRUE, "LOQ"] <- NA_real_

    #Coerce any negative Value_SD to NA
    if(!obj$settings_preprocess$suppress.messages){
      if(any((data$Value_SD<0) %in% TRUE)){
        message(paste('If value_SD < 0, replacing Value_SD with NA.',
                      sum((data$Value_SD < 0) %in% TRUE),
                      "Value_SDs will be replaced with NA.\n"))
      }
    }
    data[(data$Value_SD<0) %in% TRUE, "Value_SD"] <- NA_real_


    #Coerce any negative Time to NA
    if(!obj$settings_preprocess$suppress.messages){
      if(any((data$Time<0) %in% TRUE)){
        message(paste('If Time < 0, replacing Time with NA.',
                      sum((data$Time < 0) %in% TRUE),
                      "Times will be replaced with NA.\n"))
      }
    }
    data[(data$Time<0) %in% TRUE, "Time"] <- NA_real_

    #If any non-NA Value is currently less than its non-NA LOQ,
    #then replace it with NA
    if(!obj$settings_preprocess$suppress.messages){
      if(any((data$Value <= data$LOQ) %in% TRUE)){
        message(paste('If value <= LOQ, replacing Value with NA.',
                      sum((data$Value <= data$LOQ) %in% TRUE),
                      "Values will be replaced with NA.\n"))
      }
    }
    data$Value <- ifelse((data$Value <= data$LOQ) %in% TRUE,
                         NA_real_,
                         data$Value)

    # Impute LOQ:
    # as calc_loq_factor * minimum non-NA value in each loq_group

    if(obj$settings_preprocess$impute_loq %in% TRUE){
      if(any(is.na(data$LOQ))){
        if(!obj$settings_preprocess$suppress.messages){
          message(paste0("Estimating ",
                         sum(is.na(data$LOQ)),
                         " missing LOQs as ",
                         obj$settings_preprocess$calc_loq_factor,
                         "* minimum non-NA, non-zero Value for each unique combination of ",
                         paste(
                           sapply(obj$settings_preprocess$loq_group,
                                  as_label),
                           collapse = " + "),
                         "\n"
          )
          )
        }
        data <- do.call(dplyr::group_by,
                        args =c(list(data),
                                obj$settings_preprocess$loq_group)) %>%
          dplyr::mutate(LOQ_orig = LOQ,
                        LOQ = ifelse(is.na(LOQ_orig),
                                     {
                                       if(any((Value > 0) %in% TRUE)){
                                       min(Value[(Value > 0) %in% TRUE], na.rm = TRUE) *
                                         obj$settings_preprocess$calc_loq_factor
                                       }else{
                                         NA_real_
                                       }
                                     },
                                             LOQ_orig)
          ) %>%
          dplyr::ungroup() %>%
          as.data.frame()
      }

    } #end if impute_loq %in% TRUE

    if(!obj$settings_preprocess$suppress.messages){
      if(any((data$Value <  data$LOQ) %in% TRUE)){
      message(paste0("Converting 'Value' values of less than LOQ to NA. ",
                     sum((data$Value <  data$LOQ) %in% TRUE),
                     " values will be converted.\n"))
      }
    }
    data[(data$Value < data$LOQ) %in% TRUE, "Value"] <- NA_real_

    data$exclude <- FALSE
    data$exclude_reason <- NA_character_

    #Exclude any remaining cases where both Value and LOQ are NA
    #because if we don't have a value or an LOQ, we can't do anything
    if(any(is.na(data$Value) & is.na(data$LOQ))){
      if(!obj$settings_preprocess$suppress.messages){
        message(paste0("Excluding observations where both Value and LOQ were NA. ",
                       sum(is.na(data$Value) & is.na(data$LOQ)),
                       " observations will be excluded.\n"))
      }
      data$exclude <- ifelse(is.na(data$Value) & is.na(data$LOQ),
                             TRUE,
                             data$exclude)

      data$exclude_reason <- ifelse(is.na(data$Value) & is.na(data$LOQ),
                                    paste2(data$exclude_reason,
                                           "Value & LOQ both NA",
                                           sep = "; "),
                                    data$exclude_reason)
    }

    # For any cases where N_Subjects is NA, impute N_Subjects = 1
    if(any(is.na(data$N_Subjects))){
      data$N_Subjects_orig <- data$N_Subjects
      if(!obj$settings_preprocess$suppress.messages){
        message(
          paste0(
            "N_Subjects is NA for ",
            sum(is.na(data$N_Subjects)),
            " observations. It will be assumed = 1.\n"
          )
        )
      }

      data[is.na(data$N_Subjects), "N_Subjects"] <- 1
    }

    #for anything with N_Subjects == 1, set Value_SD to 0
    if(any(data$N_Subjects == 1)){
      if(!obj$settings_preprocess$suppress.messages){
        message(
          paste0(
            "N_Subjects is 1 for ",
            sum(data$N_Subjects == 1),
            " observations. Value_SD will be set to 0 for these observations.\n"
          )
        )
      }
      data[data$N_Subjects == 1, "Value_SD"] <- 0
    }

    # Impute missing SDs
    data$Value_SD_orig <- data$Value_SD
    if(obj$settings_preprocess$impute_sd %in% TRUE){
      if(any((data$N_Subjects >1) %in% TRUE & is.na(data$Value_SD))){
        if(!obj$settings_preprocess$suppress.messages){
          #number of SDs to be estimated
          n_sd_est <- sum(
            (data$N_Subjects >1) %in% TRUE &
              is.na(data$Value_SD)
          )
          message(paste0("Estimating missing concentration SDs (for data points with N_Subjects > 1) as ",
                         "minimum non-missing SD for each group of data given by the unique combination of variables in ",
                         paste(
                           sapply(obj$settings_preprocess$sd_group,
                                  as_label),
                           collapse = " + "),
                         ". If all SDs are missing in a group, ",
                         "SD for each observation will be imputed as",
                         # "30% of the observed mean concentration. ",
                         " 0. ",
                         n_sd_est, " missing SDs will be estimated.\n"))
        }
        data <- do.call(dplyr::group_by,
                        args = c(list(data),
                                 obj$settings_preprocess$sd_group)) %>%
          dplyr::mutate(Value_SD_orig = Value_SD,
                        Value_SD = ifelse(is.na(Value_SD_orig) &
                                                    N_Subjects > 1,
                                                  ifelse(rep(all(is.na(Value_SD_orig)), dplyr::n()),
                                                                 # 0.3 * Value,
                                                         0,
                                                                 min(Value_SD_orig, na.rm = TRUE)),
                                                  Value_SD_orig)
          ) %>%
          dplyr::ungroup() %>%
          as.data.frame()
      }

    }#end if impute_sd %in% TRUE

    #Exclude any remaining multi-subject observations where SD is NA
    if(any((data$N_Subjects >1) %in% TRUE & is.na(data$Value_SD))){
      if(!obj$settings_preprocess$suppress.messages){
        message(paste0("Excluding remaining observations with N_Subjects > 1 where reported SD is NA. ",
                       sum((data$N_Subjects >1) %in% TRUE &
                             is.na(data$Value_SD)),
                       " observations will be excluded.\n"))
      }
      data$exclude <- ifelse((data$N_Subjects >1) %in% TRUE & is.na(data$Value_SD),
                             TRUE,
                             data$exclude)
      data$exclude_reason <- ifelse((data$N_Subjects >1) %in% TRUE & is.na(data$Value_SD),
                                    paste2(data$exclude_reason,
                                            "N_Subjects > 1 and Value_SD is NA",
                                          sep = "; "),
                                    data$exclude_reason)


      # data <- subset(data,
      #                !((N_Subjects >1) %in% TRUE & is.na(Value_SD))
     # )

      # if(!obj$settings_preprocess$suppress.messages){
      #   message(paste(dim(data)[1], "observations of",
      #                 length(unique(data$Chemical)), "unique chemicals,",
      #                 length(unique(data$Species)), "unique species, and",
      #                 length(unique(data$Reference)), "unique references remain."))
      # }
    }

    #Exclude any remaining multi-subject observations where Value is NA
    if(any((data$N_Subjects >1) %in% TRUE & is.na(data$Value))){
      if(!obj$settings_preprocess$suppress.messages){
        message(paste0("Excluding observations with N_Subjects > 1 where reported Value is NA (because log-likelihood for non-detect multi-subject observations has not been implemented). ",
                       sum((data$N_Subjects >1) %in% TRUE &
                             is.na(data$Value)),
                       " observations will be excluded.\n"))
      }
      # data <- subset(data,
      #                !((N_Subjects >1) %in% TRUE & is.na(Value))
      # )

      data$exclude <- ifelse((data$N_Subjects >1) %in% TRUE & is.na(data$Value),
                             TRUE,
                             data$exclude)

      data$exclude_reason <- ifelse((data$N_Subjects >1) %in% TRUE & is.na(data$Value),
                                    paste2(data$exclude_reason,
                                            "N_Subjects > 1 and Value is NA",
                                          sep = "; "),
                                    data$exclude_reason
      )

      # if(!obj$settings_preprocess$suppress.messages){
      #   message(paste(dim(data)[1], "observations of",
      #                 length(unique(data$Chemical)), "unique chemicals,",
      #                 length(unique(data$Species)), "unique species, and",
      #                 length(unique(data$Reference)), "unique references remain.\n"))
      # }
    }

    #Exclude any NA time values
    if(any(is.na(data$Time))){
      if(!obj$settings_preprocess$suppress.messages){
        message(paste0("Excluding observations with NA time values.",
                       sum(is.na(data$Time)),
                       " observations will be excluded.\n"))
      }
      # data <- subset(data, !is.na(Time))

      data$exclude <- ifelse(is.na(data$Time),
                             TRUE,
                             data$exclude)

      data$exclude_reason <- ifelse(is.na(data$Time),
                                    paste2(data$exclude_reason,
                                            "Time is NA",
                                          sep = "; "),
                                   data$exclude_reason
      )

      # if(!obj$settings_preprocess$suppress.messages){
      #   message(paste(dim(data)[1], "observations of",
      #                 length(unique(data$Reference)), "unique references remain.\n"))
      # }
    }

    #Exclude any Dose = 0 observations
    if(any(data$Dose <= .Machine$double.eps)){
      if(!obj$settings_preprocess$suppress.messages){
        message(paste0("Excluding observations with Dose == 0 (control observations). ",
                       sum(data$Dose <= .Machine$double.eps),
                       " observations will be excluded.\n"))
      }
      # data <- subset(data, data$Dose > .Machine$double.eps)

      data$exclude <- ifelse(data$Dose < .Machine$double.eps,
                             TRUE,
                             data$exclude)

      data$exclude_reason <- ifelse(data$Dose < .Machine$double.eps,
                                    paste2(data$exclude_reason,
                                            "Dose == 0",
                                          sep = "; "),
                                   data$exclude_reason
      )

      # if(!obj$settings_preprocess$suppress.messages){
      #   message(paste(dim(data)[1], "observations of",
      #                 length(unique(data$Reference)), "unique references remain\n."))
      # }
    }

    #apply time transformation
    data$Time_orig <- data$Time
    data$Time.Units_orig <- data$Time.Units

    #first, default to identity transformation if none is specified
    if(is.null(obj$scales$time$new_units)){
      obj$scales$time$new_units <- "identity"
    }


    from_units <- unique(data$Time.Units)
    to_units <- ifelse(obj$scales$time$new_units %in% "identity",
                       from_units,
                       obj$scales$time$new_units)



    if(obj$scales$time$new_units %in% "auto"){
      to_units <- auto_units(y = data$Time,
                             from = from_units)
    }

    if(!(from_units == to_units)){
    if(!obj$settings_preprocess$suppress.messages){
      message(paste("Converting time from",
                    from_units,
                    "to",
                    to_units))
    }

    data$Time_trans <- tryCatch(convert_time(x = data$Time,
                                             from = from_units,
                                             to = to_units,
                                             inverse = FALSE),
                                error = function(err){
                                  warning(paste("invivopkfit::preprocess_data.pk():",
                                                "Error in transforming time using convert_time():",
                                                err$message))
                                  return(NA_real_)
                                })
    }else{
      data$Time_trans <- data$Time
    }
    data$Time_trans.Units <- to_units



    #Scale concentration by ratio_conc_dose

    if(is.null(obj$scales$conc$expr)){
      obj$settings_preprocess$ratio_conc_dose <- 1
    }

    ratio_conc_dose <- obj$settings_preprocess$ratio_conc_dose

    if(!obj$settings_preprocess$suppress.messages){
      if(!(ratio_conc_dose %in% 1)){
      message(paste("Concentrations, LOQs, and concentration SDs are being scaled by",
                    "ratio_conc_dose = ",
                    ratio_conc_dose, "\n"))
      }
    }


if("Conc" %in% names(data)){
      if(!obj$settings_preprocess$suppress.messages){
        message(paste("Harmonized variable `Conc` already exists in data!.",
                      "It will be overwritten as",
                      "`pmax(Value/ratio_conc_dose, LOQ/ratio_conc_dose, na.rm = TRUE)`.\n"))
      }
}
      data$Conc <- pmax(data$Value/ratio_conc_dose,
                        data$LOQ/ratio_conc_dose,
                        na.rm = TRUE)

    if("Detect" %in% names(data)){
      if(!obj$settings_preprocess$suppress.messages){
        message(paste("Harmonized variable `Detect`already exists in data!",
                      "It will be overwritten as",
                      "`is.na(Value)`.\n"))
      }
    }
      data$Detect <- !is.na(data$Value)

    if("Conc_SD" %in% names(data)){
      if(!obj$settings_preprocess$suppress.messages){
        message(paste("Harmonized variable `Conc_SD` already exists in data!",
                      "It will be overwritten as",
                      "`ifelse(data$Detect,
                        data$Value_SD/ratio_conc_dose,
                        NA_real_)\n"))
      }
    }
      data$Conc_SD <- data$Value_SD/ratio_conc_dose

      if("Conc.Units" %in% names(data)){
        if(!obj$settings_preprocess$suppress.messages){
          message(paste("Harmonized variable `Conc.Units` already exists in data!",
                        "It will be overwritten as",
                       "Value.Units/ratio_conc_dose\n"))
        }
      }

    data$Conc.Units <- ifelse(ratio_conc_dose == 1,
                              data$Value.Units,
                              paste0(data$Value.Units,
                                    "/",
                                    ratio_conc_dose))

    #apply concentration transformation.
    #
    #.conc is a placeholder that refers to any concentration variable (Conc,
    #value, Value_SD, Conc_SD, LOQ).

    #apply conc transformation
    #use tidy evaluation: Specify an expression to evaluate in the context of a data.frame
    #Transform Conc (this is either the measured value or the LOQ, depending on Detect)
    #If no conc transformation specified, assume identity
    if(is.null(obj$scales$conc$expr)){
      obj$scales$conc$dose_norm <- FALSE
      obj$scales$conc$log10_trans <- FALSE
      obj$scales$conc$expr <- rlang::new_quosure(quote(.conc),
                                                 env = caller_env())
    }

    if(!obj$settings_preprocess$suppress.messages){
      message(paste("Applying transformations to concentration variables:",
                    paste0("dose_norm = ", obj$scales$conc$dose_norm),
                    paste0("log10_trans = ", obj$scales$conc$log10_trans),
                    sep = "\n"))
      if("Conc_trans" %in% names(data)){
        message("Warning: variable `Conc_trans` already exists in the data. It will be overwritten!\n")
      }
      if("Conc_SD_trans" %in% names(data)){
        message("Warning: variable `Conc_SD_trans` already exists in the data. It will be overwritten!\n")
      }
      if("Conc_trans.Units" %in% names(data)){
        message("Warning: variable `Conc_trans.Units` already exists in the data. It will be overwritten!\n")
      }
    }

    data$Conc_trans <- rlang::eval_tidy(expr = obj$scales$conc$expr,
                                        data = cbind(data,
                                                     data.frame(.conc = data$Conc)
                                        ))
    #Transform Conc_SD
    data$Conc_SD_trans <- rlang::eval_tidy(obj$scales$conc$expr,
                                           data = cbind(data,
                                                        data.frame(.conc = data$Conc_SD)
                                           ))
    #Record new conc units
    data$Conc_trans.Units <- gsub(
      x = gsub(x = rlang::as_label(obj$scales$conc$expr),
                                  pattern = ".conc",
                                  replacement = paste0("(", unique(data$Conc.Units), ")"),
                                  fixed = TRUE),
      pattern = "Dose",
      replacement = paste0("(", unique(data$Dose.Units), ")"),
      fixed = TRUE)

    #If Series ID is not included, then assign it as NA
    if(!("Series_ID" %in% names(data))){
      data$Series_ID <- NA_integer_
    }

    # add data & data info to object
    obj$data <- data

    #if keep_data_original is FALSE, then delete original data from object, as
    #it is no longer needed
    if(obj$settings_preprocess$keep_data_original %in% FALSE){
      obj$data_original <- NULL
    }

    obj$status <- status_preprocess #preprocessing complete

    return(obj)
  }
}

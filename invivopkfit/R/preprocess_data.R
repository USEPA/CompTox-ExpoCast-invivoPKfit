#' Pre-process concentration vs. time data for analysis
#'
#' Clean data to set up for PK model fitting.
#'
#' # Preprocessing steps
#'
#' This function does the following things in the following order, and is
#' verbose about it unless told otherwise (by setting `suppress.messages =
#' FALSE`):
#'
#' - Renames columns to the standard names that invivoPKfit uses internally, and
#' adds any missing columns. See [rename_columns()].
#' - Converts concentrations to numeric, if they are not already.
#' - Converts doses to numeric, if they are not already.
#' - Converts times to numeric, if they are not already.
#' - Converts reference ID to character, if not already.
#' - Converts extraction source ID to character, if not already.
#' - If reference ID is NA, sets it to be the same as extraction source ID. (By
#' default in CvTdb, reference ID is only set if it is different from
#' extraction-source ID. For example, if data were extracted from a figure in a
#' meta-analysis paper that republished data from multiple other publications,
#' the extraction-source ID would be for the meta-analysis paper, but the
#' reference ID would be for the original publication.)
#' - Removes all observations where `DTXSID_Dosed` and `DTXSID_Analyzed` are not
#' the same (i.e., keeps only data where the parent chemical is what was
#' monitored). Creates new columns `Compound`, `DTXSID`, and `CAS` identifying
#' the chemical both dosed and analyzed.
#' - Harmonizes routes recorded as "oral" to "po" and "intravenous" to "iv".
#' - Removes all observations that have routes other than those listed in `routes_keep`.
#' - Removes all observations in media other than those listed in `media_keep`.
#' - Adds variable `iv`: a TRUE/FALSE flag, for whether route is "iv" or not.
#' - Converts `Compound` and `Species` to lower-case.
#' - Multiplies concentrations by `ratio_conc_to_dose` (*i.e.*, by the ratio
#' between the mass units for concentration and the mass units for dose).
#' - For any concentrations reported as 0, substitute NA.
#'  - Imputes LOQ for any observations missing it, using [estimate_loq()]. LOQ
#' is imputed for each unique combination of `Reference`, `DTXSID`, `Species`,
#' and `Media` as `calc_loq_factor` times the minimum *detected* `Value`.s
#'  - Substitutes NA for any concentration observations below LOQ (these are non-detects).
#'  - Removes any remaining observations where both `Value` and `LOQ` are NA.
#' For example, this situation occurs for reference/chemical/species/media
#' combinations where no LOQ was reported, and all concentrations were reported
#' as NA, so that no LOQ could be imputed.
#' - Removes any observations with NA `Time`.
#' - Adds a variable `Time.Days` where `Time` (in hours) is converted to time in days.
#'
#'
#' @param data.set A `data.frame` of concentration-time data. Preferably
#'   `pkdataset_nheerlcleaned` or `pkdataset_nheerlorig`.
#' @param model Which general model should be fit for each chemical. Presently,
#'   only "1compartment" and "2compartment" are implemented.
#' @param modelfun Either "analytic" or "full" -- whether to fit using the
#'   analytic solution to the model, or the full ODE model. Presently,
#'   "analytic" is recommended (because the analytic solution is exact and much
#'   faster).
#' @param ratio_conc_to_dose Ratio between the mass units used to report the
#'   concentration data and the mass units used to report the dose. Default 1.
#'   For example, concentration reported in ug/L and dose reported in mg/kg/day
#'   would require `ratio_conc_to_dose = 0.001`, because 1 ug/1 mg = 1e-6 g /
#'   1e-3 g = 0.001.
#' @param names_list As for [rename_columns()]: A named list where the names are
#'   the new variable names, and the values are the old variable names. If a
#'   value is NULL, then there is no old variable corresponding to that new
#'   variable. A new variable will be added, with default value as given in
#'   `defaults_list`; if no default is given in `defaults_list`, the new
#'   variable will be filled with `NA_character`.
#' @param defaults_list As for [rename_columns()]: A named list where the names
#'   are the new variable names, and the values are the default values to fill
#'   for each new variable, if the corresponding old variable name is NULL in
#'   `names_list`.
#' @param routes_keep List of routes to retain. Default `c("po", "iv")` to
#'   retain only oral and IV administration data.
#' @param media_keep List of media to retain. Default `c("blood", "plasma")` to
#'   retain only concentrations in blood and plasma.
#' @param impute_loq Logical: TRUE to impute values for missing LOQs; FALSE to
#'   leave them alone.
#' @param impute_sd Logical: TRUE to impute values for missing sample SDs for
#'   multi-subject observations; FALSE to leave them alone
#' @param suppress.messages Logical: Whether to suppress verbose messages.
#'   Default FALSE, to be verbose.
#' @return A `data.table` containing the cleaned, harmonized data, ready for
#'   model fitting. This data.table contains the following variables:
#'   \describe{
#'   \item{`Compound`}{Chemical compound name}
#'   \item{`DTXSID`}{DSSTox Substance ID,}
#'   \item{`CAS`}{Chemical CASRN}
#'   \item{`Reference`}{Study reference document ID}
#'   \item{`Species`}{Species.}
#'   \item{`Weight`}{Body weight}
#'   \item{`Weight_Units`}{Units of body weight`}
#'   \item{`Dose`}{Dose for each observation}
#'   \item{`Time`}{Time of each observation}
#'   \item{`Time.Units`}{Units of times in `Time`}
#'   \item{`Media`}{Medium in which concentrations were measured, e.g. "blood"
#'   or "plasma".}
#'   \item{`Value`}{Concentration values of each observation}
#'   \item{`Units`}{Units of concentration values in `Value`}
#'   \item{`Route`}{Route of dose administration: either "iv" (intravenous) or
#'   "po" (oral). Observations with any other routes are removed.}
#'   \item{`Extraction`}{Study extraction document ID.}
#'   \item{`LOQ`}{Limit of quantitation for concentrations}
#'   \item{`Subject`}{Individual subject identifier, if available. }
#'   \item{`iv`}{Logical TRUE/FALSE flag indicating whether `Route == "iv"`.}
#'   \item{`Time.Days`}{Time of each observation, converted to units of days.}
#'   \item{`N.Obs.Ref`}{The number of observations remaining for each unique
#'   combination of `Reference`, `DTXSID`, and `Species`, after the removal
#'   steps described in Details.}
#'   }
#'
#'   new_name = old_name

preprocess_data <- function(data.set,
                            names_list =list(
                              "Compound_Dosed" = "studies.test_substance_name_original",
                              "DTXSID_Dosed" = "chemicals_dosed.dsstox_substance_id",
                              "CAS_Dosed" = "chemicals_dosed.dsstox_casrn",
                              "Compound_Analyzed" = "series.analyte_name_original",
                              "DTXSID_Analyzed" = "chemicals_analyzed.dsstox_substance_id",
                              "CAS_Analyzed" = "chemicals_analyzed.dsstox_casrn",
                              "Reference" = "documents_reference.id",
                              "Extraction" = "documents_extraction.id",
                              "Species" = "subjects.species",
                              "Weight" ="subjects.weight_kg",
                              "Weight.Units" = NULL,
                              "Dose" = "studies.dose_level_normalized",
                              "Dose.Units" = NULL,
                              "Time" = "conc_time_values.time_hr",
                              "Time.Units" = NULL,
                              "Media" = "series.conc_medium_normalized",
                              "Value" = "conc_time_values.conc",
                              "Value.Units" = NULL,
                              "Route" = "studies.administration_route_normalized",
                              "LOQ" = "series.loq",
                              "Subject" = "subjects.id",
                              "N_Subjects" = "series.n_subjects_in_series",
                              "Value_SD" = "conc_time_values.conc_sd"),
                            defaults_list =   list(
                              "Weight.Units" = "kg",
                              "Dose.Units" = "mg/kg",
                              "Time.Units" = "h",
                              "Value.Units" = "mg/L"),
                            ratio_conc_to_dose = 1,
                            calc_loq_factor = 0.45,
                            routes_keep = c("po", "iv"),
                            media_keep = c("blood", "plasma"),
                            impute_loq = TRUE,
                            impute_sd = TRUE,
                            suppress.messages = FALSE){


  if(!suppress.messages){
    message("Renaming data columns...")
  }
  #rename columns
  data.set <- rename_columns(data.set = data.set,
                             names_list = names_list,
                             defaults_list = defaults_list)

  if(!suppress.messages){
  ### display messages describing loaded data
  message(paste(nrow(data.set), "concentration vs. time observations loaded.\n"))
  message(paste(length(unique(data.set$DTXSID_Analyzed)), "unique analytes,",
            length(unique(data.set$Species)), "unique species, and",
            length(unique(data.set$Reference)), "unique references."))
  }

  data.set$Value_orig <- data.set$Value

  ### Coerce all 'Value' values to be numeric and say so
  if (!is.numeric(data.set$Value))
  {
    value_num <- as.numeric(data.set$Value)
    old_na <- sum(is.na(data.set$Value) | !nzchar(data.set$Value))
    new_na <- sum(is.na(value_num))
    if(!suppress.messages){
      message(paste0("Column \"Value\" converted from ",
                     class(data.set$Value),
                    " to numeric. ",
                     "Pre-conversion NAs and blanks: ",
                     old_na,
                     ". Post-conversion NAs: ",
                     new_na, "."))
    }
    data.set$Value <- value_num
    rm(value_num, old_na, new_na)
  }

  ### coerce 'Dose' values to numeric and say so
  data.set$Dose_orig <- data.set$Dose
  if (!is.numeric(data.set$Dose))
  {
    dose_num <- as.numeric(data.set$Dose)
    old_na <- sum(is.na(data.set$Dose) | !nzchar(data.set$Dose))
    new_na <- sum(is.na(dose_num))
    if(!suppress.messages){
      message(paste0("Column \"Dose\" converted from ",
                     class(data.set$Dose),
      " to numeric. ",
                     "Pre-conversion NAs and blanks: ",
                     old_na,
                     ". Post-conversion NAs: ",
                     new_na, "."))
    }
    data.set$Dose <- dose_num
    rm(dose_num, old_na, new_na)
  }

  ### coerce 'Time' values to numeric and say so
  data.set$Time_orig <- data.set$Time
  if (!is.numeric(data.set$Time))
  {
    time_num <- as.numeric(data.set$Time)
    old_na <- sum(is.na(data.set$Time) | !nzchar(data.set$Time))
    new_na <- sum(is.na(time_num))
    if(!suppress.messages){
      message(paste0("Column \"Time\" converted from ",
                     class(data.set$TIme),
                     " to numeric. ",
                     "Pre-conversion NAs and blanks: ",
                     old_na,
                     ". Post-conversion NAs: ",
                     new_na, "."))
    }
    data.set$Time <- time_num
    rm(time_num, old_na, new_na)
  }


  ### Coerce all 'Reference' values to be character, and say so
  data.set$Reference_orig <- data.set$Reference
  if(!is.character(data.set$Reference)){
    if(!suppress.messages){
      message(paste0("Column \"Reference\" converted from ",
                     class(data.set$Reference),
                     " to character."))
    }
    data.set$Reference <- as.character(data.set$Reference)

  }

  ### Coerce all 'Extraction' values to be character, and say so
  if(!is.character(data.set$Extraction)){
    if(!suppress.messages){
      message(paste0("Column \"Extraction\" converted from ",
                     class(data.set$Extraction),
                     " to character."))
    }
    data.set$Extraction <- as.character(data.set$Extraction)

  }

  #If Reference is NA, set it the same as Extraction.
  data.set$Reference_orig <- data.set$Reference
  data.set[is.na(data.set$Reference),
           "Reference"] <- data.set[is.na(data.set$Reference),
                                    "Extraction"]

  #Subset to only data where the dosed and analyzed chemical are the same --
  #i.e., measuring concentration of parent chemical.
  analyte_dosed_match <- (data.set$DTXSID_Dosed == data.set$DTXSID_Analyzed) %in% TRUE
  if(any(!analyte_dosed_match)){
  if(!suppress.messages){
    message(paste("Restricting to observations where dosed DTXSID == analyzed DTXSID eliminates",
                  sum(!analyte_dosed_match),
                  "observations.")
    )
  }
  data.set <- subset(data.set,
                     analyte_dosed_match)
  }

  #Add just "DTXSID", "Compound", "CAS" columns
  if(!suppress.messages){
    message("Creating columns DTXSID, Compound, and CAS, containing the dosed & analyzed chemical identifiers.")
  }
  data.set$DTXSID <- data.set$DTXSID_Analyzed
  data.set$Compound <- data.set$Compound_Analyzed
  data.set$CAS <- data.set$CAS_Analyzed

  if(!suppress.messages){
    message(paste(dim(data.set)[1], "observations of",
                  length(unique(data.set$DTXSID)), "unique chemicals,",
                  length(unique(data.set$Species)), "unique species, and",
                  length(unique(data.set$Reference)), "unique references remain."))
  }

  # Right now code only recognizes "po" and "iv" as routes:
  ### coerce route names, 'oral' and 'intravenous', to 'po' and 'iv'
  data.set[data.set$Route %in% "oral", "Route"] <- "po"
  data.set[data.set$Route %in% "intravenous", "Route"] <- "iv"

  if(any(!(data.set$Route %in% routes_keep))){
    if(!suppress.messages){
      message(paste("Restricting to routes in",
                    paste(routes_keep, collapse = ", "),
      "eliminates",
                    sum(!(data.set$Route %in% routes_keep)), "observations."))
    }

    ### subset to data of routes in routes_keep only
    data.set <- data.set[data.set$Route %in% routes_keep, ]
    if(!suppress.messages){
      message(paste(dim(data.set)[1], "observations of",
                    length(unique(data.set$DTXSID)), "unique chemicals,",
                    length(unique(data.set$Species)), "unique species, and",
                    length(unique(data.set$Reference)), "unique references remain."))
    }
  }

  #Subset to only media in media_keep

  if(any(!(data.set$Media %in% media_keep))){
    if(!suppress.messages){
      message(paste("Restricting to media in",
                    paste(media_keep, collapse = ", "),
                    "eliminates",
                    sum(!(data.set$Media %in% media_keep)), "observations."))
    }

    ### subset to data of with media in media_keep
    data.set <- data.set[data.set$Media %in% media_keep, ]
    if(!suppress.messages){
      message(paste(dim(data.set)[1], "observations of",
                    length(unique(data.set$DTXSID)), "unique chemicals,",
                    length(unique(data.set$Species)), "unique species, and",
                    length(unique(data.set$Reference)), "unique references remain."))
    }
  }

  #set TRUE/FALSE flag for IV administration
  data.set$iv <- data.set$Route %in% "iv"

  # Harmonize the compound names:
  ### make all compound names completely lower-case
  data.set$Compound <- tolower(data.set$Compound)

  ### do the same thing for species
  data.set$Species <- tolower(data.set$Species)

  ### normalize 'Value' by ratio_conc_to_dose
  ### this makes the mass units of Value and Dose the same --
  ### e.g. mg/L and mg/kg/day
  data.set$Value <- data.set$Value * ratio_conc_to_dose
  # do the same for LOQ and SD
  data.set$LOQ <- data.set$LOQ * ratio_conc_to_dose
  data.set$Value_SD <- data.set$Value_SD * ratio_conc_to_dose

  # #Set any 0 concentrations to NA
  # if(!suppress.messages){
  #   message(paste0("Converting 'Value' values of 0 to NA.\n",
  #                  sum(data.set$Value %in% 0),
  #                  " values will be converted."))
  # }
  # data.set[data.set$Value %in% 0, "Value"] <- NA_real_

  if(impute_loq %in% TRUE){
  # Impute LOQ if it is missing
  if(any(is.na(data.set$LOQ))){
    data.set$LOQ_orig <- data.set$LOQ
    if(!suppress.messages){
      message(paste0("Estimating missing LOQs as ",
                     calc_loq_factor,
                     "* minimum detected Value for each unique combination of ",
                     "Reference, DTXSID, Media, and Species"))
    }
    data.set <- estimate_loq(dat = data.set,
                             reference_col = "Reference",
                             chem_col = "DTXSID",
                             media_col = "Media",
                             species_col = "Species",
                             value_col = "Value",
                             loq_col = "LOQ",
                             calc_loq_factor = calc_loq_factor)
  }

  if(!suppress.messages){
    message(paste0("Converting 'Value' values of less than LOQ to NA.\n",
                   sum((data.set$Value <  data.set$LOQ) %in% TRUE),
                   " values will be converted."))
  }
  data.set[(data.set$Value < data.set$LOQ) %in% TRUE, "Value"] <- NA_real_

  } #end if impute_loq %in% TRUE

  #Remove any remaining cases where both Value and LOQ are NA
  if(any(is.na(data.set$Value) & is.na(data.set$LOQ))){
    if(!suppress.messages){
      message(paste0("Removing observations where both Value and LOQ were NA.\n",
                     sum(is.na(data.set$Value) & is.na(data.set$LOQ)),
                     " observations will be removed."))
    }
    data.set <- subset(data.set,
                       !(is.na(data.set$Value) & is.na(data.set$LOQ)))

    if(!suppress.messages){
      message(paste(dim(data.set)[1], "observations of",
                    length(unique(data.set$DTXSID)), "unique chemicals,",
                    length(unique(data.set$Species)), "unique species, and",
                    length(unique(data.set$Reference)), "unique references remain."))
    }
  }

  if(impute_sd %in% TRUE){
  #Impute missing SDs
  if(any((data.set$N_Subjects >1) %in% TRUE & is.na(data.set$Value_SD))){
    data.set$Value_SD_orig <- data.set$Value_SD
    if(!suppress.messages){
      n_sd_est <- sum(
        (data.set$N_Subjects >1) %in% TRUE &
          is.na(data.set$Value_SD)
      )
      message(paste0("Estimating missing concentration SDs for multi-subject data points as ",
                     "minimum non-missing SD for each unique combination of ",
                     "Reference, DTXSID, Media, and Species.",
                     "If all SDs are missing for such a unique combination, ",
                     "SD will be imputed equal to mean. ",
                     n_sd_est, " missing SDs will be estimated."))
    }
  data.set <- estimate_conc_sd(dat = data.set,
                           reference_col = "Reference",
                           chem_col = "DTXSID",
                           media_col = "Media",
                           species_col = "Species",
                           value_col = "Value",
                           sd_col = "Value_SD",
                           n_subj_col = "N_Subjects")
  }
  }#end if impute_sd %in% TRUE

  #Remove any remaining multi-subject observations where SD is NA
  #(with imputing SD = Mean as a fallback, this will only be cases where Value was NA)
  if(any((data.set$N_Subjects >1) %in% TRUE & is.na(data.set$Value_SD))){
    if(!suppress.messages){
      message(paste0("Removing observations with N_Subjects > 1 where reported SD is NA.\n",
                     sum((data.set$N_Subjects >1) %in% TRUE &
                           is.na(data.set$Value_SD)),
                     " observations will be removed."))
    }
    data.set <- subset(data.set,
                       !((N_Subjects >1) %in% TRUE & is.na(Value_SD))
    )

    if(!suppress.messages){
      message(paste(dim(data.set)[1], "observations of",
                    length(unique(data.set$DTXSID)), "unique chemicals,",
                    length(unique(data.set$Species)), "unique species, and",
                    length(unique(data.set$Reference)), "unique references remain."))
    }
  }

  #Remove any remaining multi-subject observations where Value is NA
  if(any((data.set$N_Subjects >1) %in% TRUE & is.na(data.set$Value))){
    if(!suppress.messages){
      message(paste0("Removing observations with N_Subjects > 1 where reported Value is NA.\n",
                     sum((data.set$N_Subjects >1) %in% TRUE &
                           is.na(data.set$Value)),
                     " observations will be removed."))
    }
    data.set <- subset(data.set,
                       !((N_Subjects >1) %in% TRUE & is.na(Value))
    )

    if(!suppress.messages){
      message(paste(dim(data.set)[1], "observations of",
                    length(unique(data.set$DTXSID)), "unique chemicals,",
                    length(unique(data.set$Species)), "unique species, and",
                    length(unique(data.set$Reference)), "unique references remain."))
    }
  }


# For any cases where N_Subjects is NA, impute N_Subjects = 1
  if(any(is.na(data.set$N_Subjects))){
    data.set$N_Subjects_orig <- data.set$N_Subjects
    if(!suppress.messages){
      message(
        paste0(
          "N_Subjects is NA for ",
          sum(is.na(data.set$N_Subjects)),
          " observations. It will be assumed = 1."
        )
      )
    }

    data.set[is.na(data.set$N_Subjects), "N_Subjects"] <- 1
  }


  #convert time from hours to days

  #first remove any NA time values
  if(any(is.na(data.set$Time))){
  if(!suppress.messages){
    message(paste0("Removing observations with NA time values.\n",
            sum(is.na(data.set$Time)),
            " observations will be removed."))
  }
  data.set <- subset(data.set, !is.na(Time))

  if(!suppress.messages){
    message(paste(dim(data.set)[1], "observations of",
                  length(unique(data.set$DTXSID)), "unique chemicals,",
                  length(unique(data.set$Species)), "unique species, and",
                  length(unique(data.set$Reference)), "unique references remain."))
  }
  }

  ### convert time from hours to days
  data.set$Time.Days <- data.set$Time/24

  if(!suppress.messages){
    "Data preprocessing complete."
  }

  return(data.set)
}

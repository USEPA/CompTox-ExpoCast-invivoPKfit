#' Pre-process concentration vs. time data for analysis
#'
#' Clean data to set up for PK model fitting.
#'
#' # Preprocessing steps
#'
#' This function does the following things in the following order, and is
#' verbose about it unless told otherwise (by setting `suppress.messages = FALSE`):
#'
#' - Renames columns to the standard names that invivoPKfit uses internally, and
#' adds any missing columns. See [rename_columns()].
#' - Converts concentrations to numeric, if they are not already.
#' - Converts doses to numeric, if they are not already.
#' - Converts times to numeric, if they are not already.
#' - Converts references to character, if they are not already.
#' - Converts the data to `data.table` format.
#' - Harmonizes routes recorded as "oral" to "po" and "intravenous" to "iv".
#' - Removes all observations that have routes other than "po" or "iv".
#' - Converts all Compound and Species to lower-case.
#' - Coerces any zero concentrations to NA (nondetect).
#' - Multiplies concentrations by `ratio_conc_to_dose` (*i.e.*, by the ratio
#' between the mass units for concentration and the mass units for dose.)
#' - Adds variable `iv`: a TRUE/FALSE flag, for whether route is "iv" or not.
#' - Removes any observations with NA Time.
#' - Adds a variable where Time (in hours) is converted to time in days.
#' - Calculates the minimum number of time steps per hour for each chemical,
#' dose, and route.
#' - Removes data from any references that had fewer than 4 observations.
#' - Removes any observations where concentration was NA (nondetect) before the
#' time of peak concentration for each reference, route, DTXSID, and species.
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
#'   would require `ratio_conc_to_dose = 0.001`, because 1 ug/1 mg = 1e-6 g
#'   / 1e-3 g = 0.001.
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
#' @param LOQ_factor Numeric. Observations with concentrations less than
#'   `LOQ_factor * LOQ` will be removed. Default 2.
#' @param suppress.messages Logical: Whether to suppress verbose messages.
#'   Default FALSE, to be verbose.
#' @return A `data.table` containing the cleaned, harmonized data, ready for
#'   model fitting. This data.table contains the following variables:
#'   \describe{
#'   \item{`Compound`}{Chemical compound name, from `compound.col`}
#'   \item{`DTXSID`}{DSSTox Substance ID, from `dtxsid.col`}
#'   \item{`CAS`}{Chemical CASRN, from `cas.col`}
#'   \item{`Reference`}{Study reference ID, from `reference.col`}
#'   \item{`Species`}{Species, from `species.col` (coereced to lowercase).}
#'   \item{`Species.Weight`}{Body weight, from `species.weight.col`}
#'   \item{`Species.Weight.Units`}{Units of body weight, from
#'   `species.weight.units.col`}
#'   \item{`Dose`}{Dose, from `dose.col`, Coerced to numeric if necessary.}
#'   \item{`Time`}{Time of each observation, from `time.col`. Coerced to numeric
#'   if necessary.}
#'   \item{`Time.Units`}{Units of times in `Time`, from `time.units.col`.}
#'   \item{`Media`}{Medium in which concentrations were measured, e.g. "blood"
#'   or "plasma". From `media.col`.}
#'   \item{`Media.Units`}{"Units" of `Media`, e.g. "normalized".}
#'   \item{`Value`}{Concentration values of each observation, from `value.col`.
#'   Coerced to numeric if necessary.}
#'   \item{`Units`}{Units of concentration values in `Value`, from `units.col`.}
#'   \item{`Route`}{Route of dose administration: either "iv" (intravenous) or
#'   "po" (oral). Observations with any other routes are removed.}
#'   \item{`Source`}{Source of data, from `source.col`.}
#'   \item{`LOQ`}{Limit of quantitation for concentrations, from `loq.col`.}
#'   \item{`Subject`}{Individual subject identifier, if available. From
#'   `subject.col`.}
#'   \item{`info`}{Additional info, if any. From `info.col`.}
#'   \item{`iv`}{Logical TRUE/FALSE flag indicating whether `Route == "iv"`.}
#'   \item{`Time.Days`}{Time of each observation, converted to units of days.}
#'   \item{`Time.Steps.PerHour`}{The smallest interval between observation time
#'   points (in hours) for each unique combination of DTXSID, Dose, and Route.}
#'   \item{`N.Obs.Ref`}{The number of observations remaining for each unique
#'   value of `Reference`, after the removal steps described in Details.}
#'   }
#'
#'  new_name = old_name

preprocess_data <- function(data.set,
                            names_list =list(
                              "Compound" = "chemicals_dosed.id",
                              "DTXSID" = "chemicals_dosed.dsstox_substance_id",
                              "CAS" = "chemicals_dosed.dsstox_casrn",
                              "Reference" = "documents_reference.id",
                              "Extraction" = "documents_extraction.id",
                              "Weight" ="subjects.weight_kg",
                              "Weight.Units" = NULL,
                              "Dose" = "studies.dose_level_normalized",
                              "Dose.Units" = NULL,
                              "Time" = "conc_time_values.time_hr",
                              "Time.Units" = NULL,
                              "Media" = "series.conc_medium_normalized",
                              "Value" = "series.conc",
                              "Value.Units" = NULL,
                              "Route" = "studies.administration_route_normalized",
                              "Source" = "documents_extraction.id",
                              "LOQ" = "series.loq",
                              "Subject" = "subjects.id"),
                            defaults_list =   list(
                              "Weight.Units" = "kg",
                              "Dose.Units" = "mg/kg",
                              "Time.Units" = "h",
                              "Value.Units" = "mg/L"),
                            ratio_conc_to_dose = 1,
                            calc_loq_factor = 0.45,
                            suppress.messages = FALSE){


  if(!suppress.messages){
    message("Renaming data columns...")
  }
  #rename columns
  data.set <- rename_columns(data.set = data.set,
                             names_list = names_list,
                             defaults_list = defaults_list)

  ### number of rows in data.set
  N.PREV <- dim(data.set)[1]

  if(!suppress.messages){
  ### display messages describing loaded data
  message(paste(N.PREV, "concentration vs. time observations loaded.\n"))
  message(paste(length(unique(data.set$DTXSID)), "unique chemicals,",
            length(unique(data.set$Species)), "unique species, and",
            length(unique(data.set$Reference)), "unique references."))
  }

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
  if(!is.character(data.set$Reference)){
    data.set$Reference <- as.character(data.set$Reference)
    if(!suppress.messages){
      message(paste0("Column \"Reference\" converted from ",
                     class(data.set$Reference),
                     " to character."))
    }
  }

  # Right now code only recognizes "po" and "iv" as routes:
  ### coerce route names, 'oral' and 'intravenous', to 'po' and 'iv'
  data.set[data.set$Route %in% "oral", "Route"] <- "po"
  data.set[data.set$Route %in% "intravenous", "Route"] <- "iv"

  ### subset to data of routes 'po' and 'iv' only
  data.set <- data.set[data.set$Route %in% c("po", "iv"), ]
  if(!suppress.messages){
  message(paste("Restricting to intravenous and oral routes eliminates",
            N.PREV - dim(data.set)[1], "observations."))
  message(paste(dim(data.set)[1], "observations of",
            length(unique(data.set$DTXSID)), "unique chemicals,",
            length(unique(data.set$Species)), "unique species, and",
            length(unique(data.set$Reference)), "unique references remain."))
  }

  ### recalcuate N.PREV after removal of rows corresponding to other routes
  N.PREV <- dim(data.set)[1]

  # Harmonize the compound names:
  ### make all compound names completely lower-case
  data.set$Compound <- tolower(data.set$Compound)

  ### do the same thing for species
  data.set$Species <- tolower(data.set$Species)

  ### normalize 'Value' by ratio_conc_to_dose
  ### this makes the mass units of Value and Dose the same --
  ### e.g. mg/L and mg/kg/day
  data.set$Value <- data.set$Value * ratio_conc_to_dose

  # Impute LOQ if it is missing
  if(any(is.na(data.set$LOQ))){
    data.set <- estimate_loq(dat = data.set,
                             reference_col = "Reference",
                             chem_col = "Compound",
                             media_col = "Media",
                             species_col = "Species",
                             value_col = "Value",
                             loq_col = "LOQ",
                             calc_loq_factor = calc_loq_factor)
  }

  #Ignore data close to LOQ:
  if(!suppress.messages){
    message(paste0("Converting 'Value' values of less than ",
                   LOQ_factor,
                   " * LOQ to NA\n.",
                   sum(data.set$Value < LOQ_factor * data.set$LOQ),
                   " values will be converted."))
  }
  data.set[data.set$Value < LOQ_factor * data.set$LOQ, "Value"] <- NA_real_

  #set TRUE/FALSE flag for IV administration
  data.set$iv <- data.set$Route %in% "iv"

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

  N.PREV <- dim(data.set)[1]

  ### convert time from hours to days
  data.set$Time.Days <- data.set$Time/24


  # How many >LOQ observations do we have per chemical/species/reference?
  data_split <- split(data.set,
                      data.set[c("Reference",
                                 "DTXSID",
                                 "Species")])
  data_split2 <- lapply(data_split,
                        function(this_dat){
                          this_dat$N.Obs.Ref <- sum(!is.na(this_dat$Value))
                          return(this_dat)
                        })
  data.set <- unsplit(data_split2,
                      data.set[c("Reference",
                                 "DTXSID",
                                 "Species")])

  if(!suppress.messages){
    "Data preprocessing complete."
  }

  return(data.set)
}

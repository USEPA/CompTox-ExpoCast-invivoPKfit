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
#' adds any missing columns. See [rename_columns].
#' - Converts concentrations to numeric, if they are not already.
#' - Converts doses to numeric, if they are not already.
#' - Converts times to numeric, if they are not already.
#' - Converts references to character, if they are not already.
#' - Converts the data to `data.table` format.
#' - Harmonizes routes recorded as "oral" to "po" and "intravenous" to "iv".
#' - Removes all observations that have routes other than "po" or "iv".
#' - Converts all Compound and Species to lower-case.
#' - Coerces any zero concentrations to NA (nondetect).
#' - Multiplies concentrations by `ratio.data.to.dose` (*i.e.*, by the ratio
#' between the mass units for concentration and the mass units for dose.)
#' - Coerces any concentrations less than `LOQ_factor * LOQ` to NA (treat as
#' nondetect).
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
#' @param ratio.data.to.dose Ratio between the mass units used to report the
#'   concentration data and the mass units used to report the dose. Default 1.
#'   For example, concentration reported in ug/L and dose reported in mg/kg/day
#'   would require `ratio.data.to.dose = 0.001`, because 1 ug/1 mg = 1e-6 g
#'   / 1e-3 g = 0.001.
#' @param compound.col Column name in `data.set` that identifies chemical
#'   compound. Default "Compound".
#' @param cas.col Column name in `data.set` to identify CASRN.
#'   Default"CAS".
#' @param reference.col Column name in `data.set` to identify reference.
#'   Default "Reference".
#' @param species.col Column name in `data.set` to identify species.
#'   Default "Species".
#' @param species.default If no species column exists in `data.set`, one
#'   will be created and filled with this value.  Default NULL.
#' @param species.weight.col Column name in `data.set` to identify species
#'   weight. Default "Species.Weight".
#' @param species.weight.units.col Column name in `data.set` to identify
#'   species weigh units. Default"Species.Weight.Units".
#' @param species.weight.units.default If no species weight units column exists
#'   in `data.set`, one will be created and filled with this value. Default
#'   NULL.
#' @param dose.col Column name in `data.set` to identify dose. Default
#'   "Dose".
#' @param time.col Column name in `data.set` to identify time. Default
#'   "Time".
#' @param time.units.col Column name in `data.set` to identify time units.
#'   Default "Time.Units."
#' @param time.units.default If no time units column exists in `data.set`,
#'   one will be created and filled with this value.  Default NULL.
#' @param media.col Column name in `data.set` to identify media. Default
#'   "Media".
#' @param media.units.col Column name in `data.set` to identify media
#'   units. Default "Media.Units".
#' @param media.units.default If no media units column exists in
#'   `data.set`, one will be created and filled with this value.  Default
#'   NULL.
#' @param value.col Column name in `data.set` to identify value. Default
#'   "Value".
#' @param units.col Column name in `data.set` to identify units. Default
#'   "Units".
#' @param units.default If no units column exists in `data.set`, one will
#'   be created and filled with this value.  Default NULL.
#' @param route.col Column name in `data.set` to identify route of
#'   administration. Default "Route".
#' @param route.default If no route column exists in `data.set`, one will
#'   be created and filled with this value.  Default NULL.
#' @param source.col Column name in `data.set` to identify source. Default
#'   "Source."
#' @param source.default If no source column exists in `data.set`, one will
#'   be created and filled with this value.  Default NULL.
#' @param loq.col Column name in `data.set` to identify LOQ. Default "LOQ".
#' @param loq.default If no LOQ column exists in `data.set`, one will be
#'   created and filled with this value.  Default NULL.
#' @param subject.col Column name in `data.set` to identify subject.
#'   Default "Subject."
#' @param subject.default If no subject column exists in `data.set`, one
#'   will be created and filled with this value.  Default NULL.
#' @param info.col Column name in `data.set` to serve as info column.
#'   Default "Info".
#' @param info.default If no info column exists in `data.set`, one will be
#'   created and filled with this value.  Default NULL.
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

preprocess_data <- function(data.set,
                            ratio.data.to.dose = 1,

                            compound.col = "Compound",
                            dtxsid.col = "DTXSID",
                            cas.col = "CAS",

                            reference.col = "Reference",

                            species.col = "Species",
                            species.default = NULL,

                            species.weight.col = "Species.Weight",
                            species.weight.units.col = "Species.Weight.Units",
                            species.weight.units.default = NULL,

                            dose.col = "Dose",

                            time.col = "Time",
                            time.units.col = "Time.Units",
                            time.units.default = NULL,

                            media.col = "Media",
                            media.units.col = "Media.Units",
                            media.units.default = NULL,

                            value.col = "Value",

                            units.col = "Units",
                            units.default = NULL,

                            route.col = "Route",
                            route.default = NULL,

                            source.col = "Source",
                            source.default = NULL,

                            loq.col = "LOQ",
                            loq.default = NULL,

                            subject.col = "Subject",
                            subject.default = NULL,

                            info.col = "info",
                            info.default = NULL,

                            LOQ_factor = 2,
                            suppress.messages = FALSE){
  if(!suppress.messages){
    message("Renaming data columns...")
  }
  #rename columns to
  data.set <- rename_columns(data.set,
                             compound.col,
                             dtxsid.col,
                             cas.col,
                             reference.col,
                             species.col,
                             species.default,
                             species.weight.col,
                             species.weight.units.col,
                             species.weight.units.default,
                             dose.col,
                             time.col,
                             time.units.col,
                             time.units.default,
                             media.col,
                             media.units.col,
                             media.units.default,
                             value.col,
                             units.col,
                             units.default,
                             route.col,
                             route.default,
                             source.col,
                             source.default,
                             loq.col,
                             loq.default,
                             subject.col,
                             subject.default,
                             info.col,
                             info.default)

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

  ### coerce data.set to data.table object
  data.set <- as.data.table(data.set)

  # Right now code only recognizes "po" and "iv" as routes:
  ### coerce route names, 'oral' and 'intravenous', to 'po' and 'iv'
  data.set[Route == "oral", Route := "po"]
  data.set[Route == "intravenous", Route:="iv"]

  ### subset to data of routes 'po' and 'iv' only
  data.set <- data.set[Route %in% c("po", "iv")]
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
  data.set[, Compound := tolower(Compound)]

  ### do the same thing for species
  data.set[, Species := tolower(Species)]

  ### coerce any 'Value' values of 0 to NA

  if(!suppress.messages){
  message(paste0("Converting 'Value' values of 0 to NA.\n",
                 data.set[Value == 0, .N],
                 " values will be converted."))
  }
  data.set[Value == 0, Value := NA_real_]

  ### normalize 'Value' by ratio.data.to.dose
  ### this makes the mass units of Value and Dose the same --
  ### e.g. mg/L and mg/kg/day
  data.set[, Value := Value * ratio.data.to.dose]

  #Ignore data close to LOQ:
  if(!suppress.messages){
    message(paste0("Converting 'Value' values of less than ",
                   LOQ_factor,
                   " * LOQ to NA\n.",
                   data.set[Value < LOQ_factor * LOQ, .N],
                   " values will be converted."))
  }
  data.set[Value < LOQ_factor * LOQ, Value := NA]

  #set an iv variable to TRUE/FALSE
  data.set[Route == 'iv', iv := TRUE]
  data.set[Route !='iv', iv := FALSE]

  #convert time from hours to days

  #first remove any NA time values
  if(!suppress.messages){
    message(paste0("Removing observations with NA time values.\n",
            data.set[is.na(Time), .N],
            " observations will be removed."))
  }
  data.set <- data.set[!is.na(Time)]

  if(!suppress.messages){
  message(paste(dim(data.set)[1], "observations of",
            length(unique(data.set$DTXSID)), "unique chemicals,",
            length(unique(data.set$Species)), "unique species, and",
            length(unique(data.set$Reference)), "unique references remain."))
  }

  N.PREV <- dim(data.set)[1]

  ### convert time from hours to days
  data.set[, Time.Days := Time / 24]

  data.set[, Max.Time.Days := max(Time.Days)]

  data.set[, Time.Steps.PerHour := {
    timepts <- c(0, sort(unique(Time)))
    timediff <- diff(timepts)
    min_timediff <- min(timediff)
    1/min_timediff
  },
  by = .(DTXSID, Dose, Route)]

  # How many >LOQ observations do we have per chemical/species/reference?
  data.set[, N.Obs.Ref := sum(!is.na(Value)),
           by = .(Reference, DTXSID, Species)]
  # Not much we can do if fewer than 4 points (for instance, can't estimate Sigma):
  data.set[, Usable := N.Obs.Ref > 3,
           by = .(DTXSID, Reference, Species, Route)]
  if(!suppress.messages){
    message(paste("Restricting to references with more than three observations",
    "above LOQ eliminates",
              data.set[Usable %in% FALSE, .N],
    "observations."))
  }

  data.set <- data.set[Usable == TRUE]
if(!suppress.messages){
  message(paste(dim(data.set)[1], "observations of",
            length(unique(data.set$DTXSID)), "unique chemicals,",
            length(unique(data.set$Species)), "unique species, and",
            length(unique(data.set$Reference)), "unique references remain."))
}
  N.PREV <- dim(data.set)[1]

  #Eliminate observations for doses that are below LOQ before the peak conc.
  #Get time of peak conc for each group
  data.set[, Tpeak := Time[which.max(Value)],
           by=.(Route,Reference,DTXSID,Species)]
  data.set[, Usable := !(Time < Tpeak & is.na(Value))]

  if(!suppress.messages){
    message(paste("Eliminating observations for doses that are",
    "below LOQ before the peak conc. is reached eliminates",
              data.set[, sum(!(Usable %in% TRUE))],
    "observations."))
  }

  ### subset data to where Usable == TRUE
  data.set <- data.set[Usable %in% TRUE]

  #Remove extraneous columns
  data.set[, c("Usable",
               "Tpeak") := NULL]

  if(!suppress.messages){
  cat(paste(dim(data.set)[1], "observations of",
            length(unique(data.set$DTXSID)), "unique chemicals,",
            length(unique(data.set$Species)), "unique species, and",
            length(unique(data.set$Reference)), "unique references remain."))
  }

  if(!suppress.messages){
    "Data preprocessing complete."
  }

  return(data.set)
}

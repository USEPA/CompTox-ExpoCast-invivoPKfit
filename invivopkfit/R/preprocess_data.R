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
                            info.default = NULL){
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

  ### display messages describing loaded data
  cat(paste(N.PREV, "concentration vs. time observations loaded.\n"))
  cat(paste(length(unique(data.set$DTXSID)), "unique chemicals,",
            length(unique(data.set$Species)), "unique species, and",
            length(unique(data.set$Reference)), "unique references.\n"))

  ### Coerce all 'Value' values to be numeric and say so
  if (is.character(class(data.set$Value)))
  {
  cat("Column \"Value\" converted to numeric.\n")
  data.set$Value <- as.numeric(data.set$Value)
  }


  ### coerce 'Dose' values to numeric and say so
  if (is.character(class(data.set$Dose)))
  {
    cat("Column \"Dose\" converted to numeric.\n")
    data.set$Dose <- as.numeric(data.set$Dose)
  }

  ### coerce data.set to data.table object
  data.set <- as.data.table(data.set)

  # Right now code only recognizes "po" and "iv" as routes (my bad):
  ### coerce route names, 'oral' and 'intravenous', to 'po' and 'iv'
  data.set[Route == "oral", Route := "po"]
  data.set[Route == "intravenous", Route:="iv"]

  ### subset to data of routes 'po' and 'iv' only
  data.set <- data.set[Route %in% c("po", "iv")]
  cat(paste("Restricting to intravenous and oral routes eliminates",
            N.PREV - dim(data.set)[1], "observations.\n"))
  cat(paste(dim(data.set)[1], "observations of",
            length(unique(data.set$DTXSID)), "unique chemicals,",
            length(unique(data.set$Species)), "unique species, and",
            length(unique(data.set$Reference)), "unique references remain.\n"))

  ### recalcuate N.PREV after removal of rows corresponding to other routes
  N.PREV <- dim(data.set)[1]

  # Harmonize the compound names:
  ### make all compound names completely lower-case
  data.set[, Compound := tolower(Compound)]

  ### do the same thing for species
  data.set[, Species := tolower(Species)]

  # This way the weight units cancel (must still pay attention to denominator
  # of data to determine units for Vd):

  ### coerce any 'Value' values of 0 to NA
  data.set$Value[data.set$Value == 0] <- NA
  cat("Converting 'Value' values of 0 to NA \n")

  ### normalize 'Value' by ratio.data.to.dose
  data.set[, Value := Value * ratio.data.to.dose]

  #Ignore data close to LOQ:
  data.set[Value < 2 * LOQ, Value := NA]

  #set an iv variable to TRUE/FALSE
  data.set[Route == 'iv', iv := TRUE]
  data.set[Route !='iv', iv := FALSE]

  #convert time from hours to days

  #first remove any NA time values
  data.set <- data.set[!is.na(Time)]

  cat(paste("Requiring time to have a value != NA eliminates",
            N.PREV - dim(data.set)[1], "observations.\n"))

  cat(paste(dim(data.set)[1], "observations of",
            length(unique(data.set$DTXSID)), "unique chemicals,",
            length(unique(data.set$Species)), "unique species, and",
            length(unique(data.set$Reference)), "unique references remain.\n"))

  N.PREV <- dim(data.set)[1]

  data.set[, Time.Days := Time / 24]

  data.set[, c('Max.Time.Days',
               'Time.Steps.PerHour') := list(max(Time.Days),
                                             1 / min(
                                               diff(
                                                 c(0,
                                                   sort(
                                                     unique(Time)
                                                   )
                                                 )
                                               )
                                             )
               ),
           by = .(DTXSID, Dose, Route)]

  # How many >LOQ observations do we have per chemical/species/reference?
  data.set[, N.Obs.Ref := dim(subset(.SD, !is.na(Value)))[1], by = .(Reference, DTXSID, Species)]
  # Not much we can do if fewer than 4 points (for instance, can't estimate Sigma'):
  data.set[, Usable := N.Obs.Ref > 3, by = .(DTXSID, Reference, Species, Route)]
  data.set <- data.set[Usable == TRUE]
  cat(paste("Restricting to references with more than three observations above LOQ eliminates",
            N.PREV - dim(data.set)[1], "observations.\n"))
  cat(paste(dim(data.set)[1], "observations of",
            length(unique(data.set$DTXSID)), "unique chemicals,",
            length(unique(data.set$Species)), "unique species, and",
            length(unique(data.set$Reference)), "unique references remain.\n"))
  N.PREV <- dim(data.set)[1]

  data.set[,
           Usable := ifelse(Time >= .SD[Value == max(Value, na.rm = TRUE), Time] |
                              !is.na(Value), TRUE, FALSE),
           by=.(Route,Reference,DTXSID,Species)]

  ### subset data to where Usable == TRUE
  data.set <- data.set[Usable == TRUE]
  cat(paste("Eliminating observations for doses that are below LOQ before the peak conc. is reached eliminates",
            N.PREV - dim(data.set)[1], "observations.\n"))
  cat(paste(dim(data.set)[1], "observations of",
            length(unique(data.set$DTXSID)), "unique chemicals,",
            length(unique(data.set$Species)), "unique species, and",
            length(unique(data.set$Reference)), "unique references remain.\n"))
}

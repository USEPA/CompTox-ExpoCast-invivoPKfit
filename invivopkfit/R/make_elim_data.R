#' Get elimination-phase data
#'
#' A helper function used by `get_starts` to get elimination-phase data to
#' roughly estimate elimination rate, volume of distribution, etc. based on
#' linear fits to elimination-phase data
#'
#' Removes any rows with `is.na(Value)` and/or `Dose == 0`. If IV data
#' available, keeps only that data. If IV data are not available but oral data are
#' available, keeps only the oral data.
#'
#' @param fitdata A concentration-time-dose data.frame
#' @return A data.frame with the subset of `fitdata` that reflects elimination.
#'   Includes new variable `ValueDose`, which is `Value` divided by `Dose`.

make_elim_data <- function(fitdata){
  if ("iv" %in% fitdata$Route){
    #For a moment, assume 1-compartment model
    #Conc/Dose = (1/Vdist)*exp(-kel*t)
    #log(Conc/Dose) = log(1/Vdist) + -kel*t
    #of form y = b + m*x, so linear regression works to fit slope of -kel

    #take IV data with valid concentrations
    #exclude dose = 0 because we will be dividing by dose
    elim_data <- subset(fitdata, Route == "iv" &
                          !is.na(Value) &
                          Dose > 0)
    #divide concentration by dose to get LHS of equation
    elim_data$ValueDose <- elim_data$Value / elim_data$Dose

  }else if ("po" %in% fitdata$Route) {
    #if we don't have IV data but we do have oral data,
    #then fit kelim from that
    elim_data <- subset(fitdata, Route == "po" &
                          !is.na(Value) &
                          Dose > 0)
    #again, divide conc by dose to get LHS of 1-compartment model
    elim_data$ValueDose <- elim_data$Value/elim_data$Dose
    #for oral data, consider only the data after absorption is complete
    #i.e. after peak plasma concentration has been reached
    max_time <- elim_data[order(-elim_data$Value, #sort conc descending
                                elim_data$Time), #then sort time ascending within conc
                          "Time"][1]

    #take only data after the time of Cmax
    elim_data <- subset(elim_data, Time >= max_time)
  } #end if ("iv" %in% fitdata$Route)/else if ("po" %in% fitdata$Route)

  return(elim_data)
}

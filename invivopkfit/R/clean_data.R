clean_data <- function(data, ratio_data_to_dose) {

  ### make conc numeric
  if (!is.numeric(class(data$conc))) {

    cat("Column \"conc\" converted to numeric. \n")

    data[, conc := as.numeric(conc)]
  }


  ### get number of rows in data
  n_prev <- dim(data)[1]

  ### display messages describing loaded data
  cat(paste(n_prev, "concentration vs time observations loaded.\n"))
  cat(paste(length(unique(data$dsstox_casrn)), "unique chemicals,",
            length(unique(data$species)), "unique species, and",
            length(unique(data$reference)), "unique references.\n"))

  ### make dose numeric
  if (!is.numeric(class(data$dose))) {

    cat("Column \"dose\" converted to numeric. \n")

    data[, dose := as.numeric(dose)]
  }

  ### change route names
  data[route == "oral", route := "po"]
  data[route == "intravenous", route := "iv"]

  ### harmonize the compound names
  data[, compound := tolower(compound)]

  ### subset to data of routes 'po' and 'iv' only
  data <- data[route %in% c("po", "iv")]

  cat(paste("Restricting to intravenous and oral routes eliminates",
            n_prev - dim(data)[1], "observations.\n"))
  cat(paste(dim(data)[1], "observations of",
            length(unique(data$cas)), "unique chemicals,",
            length(unique(data$species)), "unique species, and",
            length(unique(data$reference)), "unique references remain.\n"))

  ### harmonize the species names
  data[, species := tolower(species)]

  ### change concentrations of 0 to NA
  # because 0 should mean non-detect
  if(any(data$conc) == 0) {
  data$conc[data$conc == 0] <- NA
  cat("Converting concentrations of \"0\" to NA")
  }

  ### normalize concentrations by ratio_data_to_dose
  data[, conc := conc * ratio_data_to_dose]

  ### ignore data close to LOQ
  data[conc < 2 * loq, conc := NA]

  ### set an iv variable to TRUE/FALSE
  data[route == "iv", iv := TRUE]
  data[route != "iv", iv := FALSE]

  ### convert time from hours to days
  data <- data[!is.na(time)]
  cat(paste("Requiring time to have a value != NA eliminates",
            n_prev - dim(data)[1], "observations.\n"))
  cat(paste(dim(data)[1], "observations of",
            length(unique(data$cas)), "unique chemicals,",
            length(unique(data$species)), "unique species, and",
            length(unique(data$reference)), "unique references remain.\n"))
  n_prev <- dim(data)[1]
  data[, time_days := time / 24]
  data[, c('max_time_days',
               'time_steps_per_hour') := list(max(time_days),
                                             1 / min(diff(c(0, sort(unique(time)))))),
           by = .(cas, dose, route)]

  ### how many >LOQ observations do we have per chemical/species/reference?
  data[, n_obs_ref := dim(subset(.SD, !is.na(conc)))[1], by = .(reference, cas, species)]

  ### not much we can do if fewer than 4 points (for instance, can't estimate sigma'):
  data[, usable := n_obs_ref > 3, by = .(cas, reference, species, route)]

  ### subset to usable data
  data <- data[usable == TRUE]

  cat(paste("Restricting to references with more than three observations above LOQ eliminates",
            n_prev - dim(data)[1], "observations.\n"))
  cat(paste(dim(data)[1], "observations of",
            length(unique(data$cas)), "unique chemicals,",
            length(unique(data$species)), "unique species, and",
            length(unique(data$reference)), "unique references remain.\n"))
  n_prev <- dim(data)[1]

  ###
  data[, usable := ifelse(time >= .SD[conc == max(conc, na.rm = TRUE), time] | !is.na(conc),
                          TRUE, FALSE),
       by = .(route, reference, cas, species)]

  ### subset data to where Usable == TRUE
  data <- data[usable == TRUE]
  cat(paste("Eliminating observations for doses that are below LOQ before the peak conc. is reached eliminates",
            n_prev - dim(data)[1], "observations.\n"))
  cat(paste(dim(data)[1], "observations of",
            length(unique(data$cas)), "unique chemicals,",
            length(unique(data$species)), "unique species, and",
            length(unique(data$reference)), "unique references remain.\n"))

  n_prev <- dim(data)[1]

  return(data)
}

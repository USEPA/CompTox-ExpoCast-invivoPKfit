#' Get non-compartmental analysis (NCA) results
#'
#' Do NCA on one set of dose-normalized concentration vs. time data
#'
#' This function calls [PK::nca()] after automatically determining the correct
#' "design" expected by [PK::nca()] (depending on whether all subjects have
#' measurements at all time points or not).
#'
#' @param obs_data A `data.table`: a set of concentration vs. time data (one
#'   chemical, species, route of dose administration, and medium of
#'   concentration). Must have variables named `Time` (times of observations),
#'   `Subject` (subject ID for each observation), `Conc_Dose` (dose-normalized
#'   concentration, with non-detects substituted by LOQ), and `Route` ("po" for
#'   oral administration, "iv" for IV administration).
#' @return A named list of NCA-estimated parameters as from [PK::nca()]. These
#'   will be "AUC_tlast_1mgkg" (AUC to the last time point), "AUC_inf_1mgkg"
#'   (AUC extrapolated to infinity), "AUMC_inf_1mgkg" (AUMC to infinity), "MRT"
#'   (mean residence time), "halflife" (non-compartmental halfife), "Clearance",
#'   (clearance), and "Vss" (volume of distribution at steady-state). Note that
#'   if route is oral, halflife and Vss are not valid estimates; Clearance is
#'   actually Clearance divided by bioavailability (Fgutabs); and MRT is
#'   actually MTT.
#'
do_nca <- function(obs_data,
                    dose_norm = TRUE){
  obs_data <- as.data.frame(obs_data)
  if(all(c("Subject",
           "Study_ID",
           "Series_ID")
         %in% names(obs_data))){
  #Ensure subject IDs go with study IDs
  obs_data$Subject_ID <- paste(obs_data$Subject,
                               obs_data$Study_ID,
                               obs_data$Series_ID)
  }else{ #if no subject ID info:
    #assume only one subject
    obs_data$Subject_ID <- "1"
  }
  #First: determine whether design is "complete", "ssd", or "batch"
  #if all time points are available for all subjects, it's "complete"
  #if one measurement per subject, it's "ssd" (unless there is only one subject, in which case it's "complete")
  #if multiple time points are measured for each subject, it's "batch"
  #create a table of number of measurements per subject and time
  ntab <- with(obs_data, table(Time, Subject_ID))
  m <- matrix(ntab, ncol = ncol(ntab))
  if(ncol(m)==1){
    design <- "complete"
  }else{
    if(all(m[!diag(nrow(m))] == 0)){
      #one measurement per subject per time
      design <- "ssd"
    }else{
      if(all(m==1)){
        design <- "complete"
      }else{
        design <- "batch"
      }
    }
  }

  if(!("Conc_Dose" %in% names(obs_data))){
    obs_data$Conc_Dose <- obs_data$Conc/obs_data$Dose
  }

  #create a data frame
  if(dose_norm %in% TRUE){
    #with dose-normalized concentrations, time, and subject ID
    nca_dat <- obs_data[c("Conc_Dose", "Time", "Subject_ID")]
    nca_dose <- 1
  }else{
    nca_dat <- obs_data[c("Conc", "Time", "Subject_ID")]
    nca_dose <- unique(obs_dat$Dose)
  }
  names(nca_dat) <- c("conc", "time", "id")
  nca_dat$group <- 1

  nca_est <- tryCatch(
    suppressMessages(PK::estimator(PK::nca(data = nca_dat,
                                           dose = nca_dose,
                                           design= design,
                                           method = "z"))[, 1]),
    error = function(err) rep(NA_real_, 7)
  )


    names(nca_est) <- c("AUC_tlast",
                        "AUC_inf",
                        "AUMC_inf",
                        "MRT",
                        "halflife",
                        "Clearance",
                        "Vss")

  return(as.data.frame(as.list(nca_est)))
}

get_auc <- function(obs_data){
  #first interpolate
  interp_dat <- approx(x = obs_data$Time,
                       y = obs_data$Conc_Dose,
                       xout = sort(unique(obs_data$Time)),
                       method = "linear",
                       rule = 2,
                       ties = mean)
  if(!(any(interp_dat$x %in% 0))){
    x <- c(0,interp_dat$x)
    y <- c(0, interp_dat$y)
  }else{
    x <- interp_dat$x
    y <- interp_dat$y
  }

  if(length(x) != length(unique(x))) browser()
  #then use trapezoidal rule
  auc_tlast <-  pracma::trapz(x = x,
                              y = y)

  return(auc_tlast)
}

#' Get summary data info and NCA
#'
#' @param obj An object of class [pk()]
#' @export
#' @author Caroline Ring
data_info.pk <- function(obj){
  #check status
  objname <- deparse(substitute(obj))
  status <- obj$status
  if(status >= status_data_info){
    warning(paste0(objname,
                   " current status is ",
                   status,
                   ". data_info() will reset its status to ",
                   status_data_info,
                   ". Any results from later workflow stages will be lost."))
  }

  #if preprocessing not already done, do it
  if(obj$status < status_preprocess){
    obj <- preprocess_data(obj)
  }

  data <- obj$data
  #get the summary data info
  #unique Chemical, Species, References, Studies, Routes
  dat_info <- as.list(unique(data[c("Chemical",
                                    "Species")]))

  dat_info$References_Analyzed <- sort(unique(data$Reference))
  #get a list of studies analyzed
  dat_info$Studies_Analyzed <- sort(unique(data$Study))
  #get a list of routes analyzed
  dat_info$Routes_Analyzed <- sort(unique(data$Route))

  #get a list of media analyzed
  dat_info$Media_Analyzed <- sort(unique(data$Media))

  #get the number of detects and non-detects by route and medium
  dat_info$n_dat <- aggregate(x = list(Detect = data$Detect),
                              by = data[c("Route", "Media")],
                              FUN = function(Detect){
                                c("Detect" = sum(Detect %in% TRUE),
                                  "NonDetect" = sum(Detect %in% FALSE))
                              })
  names(dat_info$n_dat) <- gsub(x = names(dat_info$n_dat),
                                pattern = "^Detect",
                                replacement = "")

  #get time of last detected observation
  if(any(data$Detect %in% TRUE)){
    dat_info$last_detect_time <- max(data[data$Detect %in% TRUE, "Time_trans"])
  }else{
    dat_info$last_detect_time <- 0
  }
  attr(dat_info$last_detect_time, "units") <- unique(data$Time_trans.Units)

  #get time of last observation
  dat_info$last_time <- max(data$Time_trans)
  attr(dat_info$last_time, "units") <- unique(data$Time_trans.Units)

  #do NCA
  dat_info$nca <- do.call(dplyr::group_by,
                          args =c(list(data),
                                  obj$settings_data_info$nca_group)) %>%
    dplyr::summarise(tlast = max(Time_trans),
                     calc_nca(time = Time_trans,
                              dose = Dose,
                              conc = Conc,
                              detect = Detect,
                              n_subj = N_Subjects,
                              subject_id = Subject_ID)) %>%
    as.data.frame()

  # #add units for each NCA param
  dat_info$nca[dat_info$nca$param_name %in% c("AUC to tlast",
                                              "AUC to infinity"),
               "param_units"] <- paste(unique(data$Conc.Units),
                                            "*",
                                            unique(data$Time_trans.Units))
  dat_info$nca[dat_info$nca$param_name %in% "AUMC to infinity",
               "param_units"] <- paste(unique(data$Conc.Units),
                                            "*",
                                            unique(data$Time_trans.Units),
                                          "*",
                                          unique(data$Time_trans.Units))
  dat_info$nca[dat_info$nca$param_name %in% c("Mean residence time",
                                              "non-compartmental half-life",
                                              "tmax"),
               "param_units"] <-  unique(data$Time_trans.Units)
  dat_info$nca[dat_info$nca$param_name %in% "Clearance",
               "param_units"] <- paste0("1/",
                                           unique(data$Time_trans.Units))
  dat_info$nca[dat_info$nca$param_name %in% "Volume of distribution at steady state",
               "param_units"] <- paste0(unique(data$Conc.Units),
                                           "/",
                                           unique(data$Dose.Units))
  dat_info$nca[dat_info$nca$param_name %in% "Cmax",
               "param_units"] <- unique(data$Conc.Units)
  #save data summary info
  obj$data_info <- dat_info

  obj$status <- status_data_info #data summarization complete

  return(obj)
}

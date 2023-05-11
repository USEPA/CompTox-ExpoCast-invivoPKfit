data_info.pk <- function(obj){
  #check status
  objname <- deparse(substitute(obj))
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
    dat_info$last_detect_time <- max(data[data$Detect %in% TRUE, "Time"])
  }else{
    dat_info$last_detect_time <- 0
  }

  #get time of last observation
  dat_info$last_time <- max(data$Time)

  #save data summary info
  obj$data_info <- dat_info

  #do NCA

  obj$status <- 3 #data summarization complete
}

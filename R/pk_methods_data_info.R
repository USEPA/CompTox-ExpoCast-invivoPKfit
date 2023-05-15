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

  data_summary <- do.call(dplyr::group_by,
          args =c(list(data),
                  obj$settings_data_info$nca_group)) %>%
    dplyr::summarise(group_id = dplyr::cur_group_id(),
                     group_n_obs = dplyr::n(),
                     group_n_detect = sum(Detect %in% TRUE),
                     tlast = max(Time_trans),
                     tlast_detect = max(Time_trans[Detect %in% TRUE]),
                     Time.Units = unique(Time.Units),
                     Time_trans.Units = unique(Time_trans.Units),
                     Conc.Units = unique(Conc.Units),
                     Conc_trans.Units = unique(Conc_trans.Units),
                     Dose.Units = unique(Dose.Units)) %>%
    as.data.frame()

  dat_info <- list("data_summary" = data_summary)

  #do NCA
  dat_info$nca <- do.call(dplyr::group_by,
                          args =c(list(data),
                                  obj$settings_data_info$nca_group)) %>%
    dplyr::summarise(group_id = dplyr::cur_group_id(),
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

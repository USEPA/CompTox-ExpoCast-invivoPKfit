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
    dplyr::summarise(nca_group_id = dplyr::cur_group_id(),
                     nca_group_n_obs = dplyr::n(),
                     nca_group_n_detect = sum(Detect %in% TRUE),
                     tlast = max(Time_trans),
                     tlast_detect = max(Time_trans[Detect %in% TRUE]),
                     tfirst = min(Time_trans),
                     tfirst_detect = min(Time_trans[Detect %in% TRUE]),
                     Time.Units = unique(Time.Units),
                     Time_trans.Units = unique(Time_trans.Units),
                     Conc.Units = unique(Conc.Units),
                     Conc_trans.Units = unique(Conc_trans.Units),
                     Dose.Units = unique(Dose.Units)) %>%
    as.data.frame()

  dat_info <- list("data_summary" = data_summary)

  #do NCA -- only on detects.
  dat_info$nca <- do.call(dplyr::group_by,
                          args =c(list(data),
                                  obj$settings_data_info$nca_group)) %>%
    dplyr::summarise(nca_group_id = dplyr::cur_group_id(),
                     calc_nca(time = Time_trans[Detect %in% TRUE],
                              dose = Dose[Detect %in% TRUE],
                              conc = Conc[Detect %in% TRUE],
                              detect = Detect[Detect %in% TRUE],
                              route = unique(Route),
                              n_subj = N_Subjects[Detect %in% TRUE],
                              subject_id = Subject_ID[Detect %in% TRUE])) %>%
    as.data.frame()

  # #add units for each NCA param
  dat_info$nca[dat_info$nca$param_name %in% c("AUC_tlast",
                                              "AUC_infinity"),
               "param_units"] <- paste(unique(data$Conc.Units),
                                            "*",
                                            unique(data$Time_trans.Units))
  dat_info$nca[dat_info$nca$param_name %in% "AUMC_infinity",
               "param_units"] <- paste(unique(data$Conc.Units),
                                            "*",
                                            unique(data$Time_trans.Units),
                                          "*",
                                          unique(data$Time_trans.Units))
  dat_info$nca[dat_info$nca$param_name %in% c("MRT",
                                              "MTT",
                                              "halflife",
                                              "tmax"),
               "param_units"] <-  unique(data$Time_trans.Units)
  dat_info$nca[dat_info$nca$param_name %in% c("CLtot",
                                              "CLtot/Fgutabs"),
               "param_units"] <- paste0("1/",
                                           unique(data$Time_trans.Units))
  dat_info$nca[dat_info$nca$param_name %in% "Vss",
               "param_units"] <- paste0(unique(data$Conc.Units),
                                           "/",
                                           unique(data$Dose.Units))
  dat_info$nca[dat_info$nca$param_name %in% "Cmax",
               "param_units"] <- unique(data$Conc.Units)
  #save data summary info
  obj$data_info <- dat_info

  #get data flags:
  data_flags <- NULL
  df <- merge(dat_info$data_summary,
              dat_info$nca,
              by = intersect(names(dat_info$data_summary),
                             names(dat_info$nca))
  )
  #for groups with route == "oral", is Cmax equal to the first or last detected time?
  tmax_first <- df %>%
    dplyr::filter(Route %in% "oral" &
                    param_name %in% "tmax") %>%
                    dplyr::mutate(test = isTRUE(all.equal(param_value, tfirst_detect,
                                            tolerance = sqrt(.Machine$double.eps))))

  if(sum(tmax_first$test)>0){
      data_flags <- c(data_flags,
                      paste("tmax is equal to time of first detect in",
                            sum(tmax_first$test),
                            "NCA data groups with oral data"))
  }

  tmax_last <- df %>%
    dplyr::filter(Route %in% "oral" &
                    param_name %in% "tmax") %>%
    dplyr::mutate(test = isTRUE(all.equal(param_value, tlast_detect,
                                     tolerance = sqrt(.Machine$double.eps))))

  if(sum(tmax_last$test)>0){
    data_flags <- c(data_flags,
                    paste("tmax is equal to time of last detect in",
                          sum(tmax_last$test),
                          "NCA data groups with oral data"))
  }

  #Is the NCA clearance negative?

  neg_CLtot <- df %>%
    dplyr::filter(param_name %in% c("CLtot",
                                    "CLtot/Fgutabs") &
                    param_value < 0)

  if(nrow(neg_CLtot)>0){
    data_flags <- c(data_flags,
                    paste("CLtot or CLtot/Fgutabs is negative in",
                          nrow(neg_CLtot),
                          "NCA data groups"))
  }

  neg_AUCinf <- df %>%
    dplyr::filter(param_name %in% c("AUC_infinity") &
                    param_value < 0)

  if(nrow(neg_AUCinf)>0){
    data_flags <- c(data_flags,
                    paste("AUC_infinity is negative in",
                          nrow(neg_AUCinf),
                          "NCA data groups"))
  }

  obj$data_info$data_flags <- data_flags

  obj$status <- status_data_info #data summarization complete

  return(obj)
}

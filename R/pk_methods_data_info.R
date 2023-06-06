#' Get summary data info and NCA
#'
#' @param obj An object of class [pk()]
#' @export
#' @importFrom magrittr `%>%`
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
                   ". Any results from later workflow stages will be lost.\n"))
  }

  #if preprocessing not already done, do it
  if(obj$status < status_preprocess){
    obj <- preprocess_data(obj)
  }

  if(obj$settings_preprocess$suppress.messages %in% FALSE){
    message("data_info.pk(): Getting data summary info\n")
  }

  data <- obj$data

  data_summary <- do.call(dplyr::group_by,
          args =c(list(data),
                  obj$settings_data_info$nca_group)) %>%
    dplyr::summarise(nca_group_id = dplyr::cur_group_id(),
                     n_obs = dplyr::n(),
                     n_exclude = sum(exclude %in% TRUE),
                     n_detect = sum(Detect %in% TRUE & exclude %in% FALSE),
                     n_series_id = length(unique(Series_ID[exclude %in% FALSE])),
                     n_timepts = length(unique(Time_trans[exclude %in% FALSE])),
                     tlast = ifelse(any(exclude %in% FALSE),
                                    max(Time_trans[exclude %in% FALSE]),
                                    NA_real_),
                     tlast_detect = ifelse(any(Detect %in% TRUE & exclude %in% FALSE),
                                           max(Time_trans[Detect %in% TRUE & exclude %in% FALSE]),
                                           NA_real_),
                     tfirst = ifelse(any(exclude %in% FALSE),
                                     min(Time_trans[exclude %in% FALSE]),
                                     NA_real_),
                     tfirst_detect = ifelse(any(Detect %in% TRUE & exclude %in% FALSE),
                                            min(Time_trans[Detect %in% TRUE & exclude %in% FALSE]),
                                            NA_real_),
                     Time.Units = unique(Time.Units),
                     Time_trans.Units = unique(Time_trans.Units),
                     Conc.Units = unique(Conc.Units),
                     Conc_trans.Units = unique(Conc_trans.Units),
                     Dose.Units = unique(Dose.Units)) %>%
    as.data.frame()



  #do NCA
  if(obj$settings_preprocess$suppress.messages %in% FALSE){
    message("data_info.pk(): Doing non-compartmental analysis\n")
  }

  nca <- do.call(dplyr::group_by,
                          args =c(list(data),
                                  obj$settings_data_info$nca_group)) %>%
    dplyr::summarise(nca_group_id = dplyr::cur_group_id(),
                     Conc.Units = unique(Conc.Units),
                    Time_trans.Units = unique(Time_trans.Units),
                    Dose.Units = unique(Dose.Units),
                     calc_nca(time = Time_trans[exclude %in% FALSE],
                              dose = Dose[exclude %in% FALSE],
                              conc = Conc[exclude %in% FALSE],
                              detect = Detect[exclude %in% FALSE],
                              route = unique(Route),
                              series_id = Series_ID[exclude %in% FALSE])) %>%
    dplyr::group_by(nca_group_id) %>%
    dplyr::mutate(param_units = dplyr::case_when(
      param_name %in% c("AUC_tlast",
                        "AUC_infinity") ~ paste(Conc.Units,
                                                "*",
                                                Time_trans.Units),
      param_name %in% "AUMC_infinity" ~ paste(Conc.Units,
                                              "*",
                                              Time_trans.Units,
                                              "*",
                                              Time_trans.Units),
      param_name %in% c("MRT",
                        "MTT",
                        "halflife",
                        "tmax") ~ Time_trans.Units,
      param_name %in% c("CLtot",
                        "CLtot/Fgutabs") ~ paste0("L/",
                                                  Time_trans.Units),
      param_name %in% "Vss" ~ paste0(Conc.Units,
                                     "/",
                                     Dose.Units),
      param_name %in% "Cmax" ~ Conc.Units
      )) %>%
    dplyr::select(-c(Conc.Units,
                     Time_trans.Units,
                     Dose.Units)) %>%
    as.data.frame()

  #get data flags:
  if(obj$settings_preprocess$suppress.messages %in% FALSE){
    message("data_info.pk(): Getting data flags\n")
  }

  nca_wide <- nca %>% tidyr::pivot_wider(id_cols = nca_group_id,
                                 names_from = param_name,
                                 values_from = param_value)

  df <- dplyr::inner_join(data_summary,
              nca_wide,
              by = "nca_group_id") %>%
    dplyr::mutate(`Cmax/Dose` = Cmax/Dose,
                  `AUC/Dose` = AUC_infinity/Dose) %>%
    dplyr::mutate(
    data_flag = ifelse(
      n_obs - n_exclude < 3,
      "Fewer than 3 non-excluded observations",
      NA_character_
    )
  ) %>% dplyr::mutate(
    data_flag = ifelse(
      n_detect < 3,
      paste2(data_flag,
             "Fewer than 3 detected observations",
             sep = " | "),
      data_flag
    )
  ) %>% dplyr::mutate(
    data_flag = ifelse(
      Route %in% "oral" &
        isTRUE(all.equal(tmax, tfirst_detect,
                         tolerance = sqrt(.Machine$double.eps))) &
                 !is.na(tmax),
      paste2(data_flag,
             "tmax is equal to time of first detect",
             sep = " | "),
      data_flag
    )) %>% dplyr::mutate(
      data_flag = ifelse(
        Route %in% "oral" &
          isTRUE(all.equal(tmax, tlast_detect,
                           tolerance = sqrt(.Machine$double.eps))) &
                   !is.na(tmax),
        paste2(data_flag,
               "tmax is equal to time of last detect",
               sep = " | "),
        data_flag
      )) %>% dplyr::mutate(
        data_flag = ifelse(
          (CLtot < 0) %in% TRUE |
            (`CLtot/Fgutabs` < 0) %in% TRUE,
          paste2(data_flag,
                 "CLtot or CLtot/Fgutabs is negative",
                 sep = " | "),
          data_flag)
      ) %>% dplyr::mutate(
        data_flag = ifelse((AUC_infinity < 0) %in% TRUE,
                           paste2(data_flag,
                                  "AUC_infinity is negative",
                                  sep = " | "),
                           data_flag)
      ) %>%
    as.data.frame()

  #do NCA on all the data together, grouped by Chemical, Species, Route, Media
  nca_dose_norm <- data %>% dplyr::group_by(Chemical, Species, Route, Media) %>%
    dplyr::summarise(nca_dosenorm_group_id = dplyr::cur_group_id(),
                     n_obs = dplyr::n(),
                     n_exclude = sum(exclude %in% TRUE),
                     n_detect = sum(Detect %in% TRUE & exclude %in% FALSE),
                     n_series_id = length(unique(Series_ID[exclude %in% FALSE])),
                     n_timepts = length(unique(Time_trans[exclude %in% FALSE])),
                     n_ref = length(unique(Reference[exclude %in% FALSE])),
                     tlast = ifelse(any(exclude %in% FALSE),
                                    max(Time_trans[exclude %in% FALSE]),
                                    NA_real_),
                     tlast_detect = ifelse(any(exclude %in% FALSE),
                                           max(Time_trans[Detect %in% TRUE & exclude %in% FALSE]),
                                           NA_real_),
                     tfirst = ifelse(any(exclude %in% FALSE),
                                     min(Time_trans[exclude %in% FALSE]),
                                     NA_real_),
                     tfirst_detect = ifelse(any(exclude %in% FALSE),
                                            min(Time_trans[Detect %in% TRUE & exclude %in% FALSE]),
                                            NA_real_),
                       Conc.Units = unique(Conc.Units),
                     Time_trans.Units = unique(Time_trans.Units),
                     Dose.Units = unique(Dose.Units),
                     calc_nca(time = Time_trans[exclude %in% FALSE],
                              dose = 1,
                              conc = Conc[exclude %in% FALSE]/Dose[exclude %in% FALSE],
                              detect = Detect[exclude %in% FALSE],
                              route = unique(Route),
                              series_id = Series_ID[exclude %in% FALSE])) %>%
    dplyr::mutate(param_units = dplyr::case_when(
      param_name %in% c("AUC_tlast",
                        "AUC_infinity") ~ paste(Conc.Units,
                                                "*",
                                                Time_trans.Units),
      param_name %in% "AUMC_infinity" ~ paste(Conc.Units,
                                              "*",
                                              Time_trans.Units,
                                              "*",
                                              Time_trans.Units),
      param_name %in% c("MRT",
                        "MTT",
                        "halflife",
                        "tmax") ~ Time_trans.Units,
      param_name %in% c("CLtot",
                        "CLtot/Fgutabs") ~ paste0("L/",
                                                  Time_trans.Units),
      param_name %in% "Vss" ~ paste0(Conc.Units,
                                     "/",
                                     Dose.Units),
      param_name %in% "Cmax" ~ Conc.Units
    )) %>%
    dplyr::select(-c(Conc.Units,
                     Dose.Units,
                     Time_trans.Units)) %>%
    as.data.frame()

  #Other data flags:
  #Check for obeying dose normalization by route

  dose_norm_check <- df %>%
    dplyr::group_by(Chemical, Species, Route, Media) %>%
    dplyr::summarise(
      Cmax_Dose_fold_range = {
        tmprange <- suppressWarnings(range(`Cmax/Dose`, na.rm = TRUE))
        if(all(!is.finite(tmprange))) tmprange <- c(NA_real_, NA_real_)
        tmprange[2]/tmprange[1]
      },
      data_flag_Cmax = ifelse(Cmax_Dose_fold_range > 2,
                         "Cmax may not scale with dose. Cmax/Dose range > 2-fold across NCA groups for this Route/Media",
                         NA_character_),
      AUC_Dose_fold_range = {
        tmprange <- suppressWarnings(range(`AUC/Dose`, na.rm = TRUE))
        if(all(!is.finite(tmprange))) tmprange <- c(NA_real_, NA_real_)
        tmprange[2]/tmprange[1]
      },
      data_flag_AUC = ifelse(AUC_Dose_fold_range > 2,
                              "AUC_infinity may not scale with dose. AUC/Dose range > 2-fold across NCA groups for this Route/Media",
                              NA_character_),
    ) %>% as.data.frame()


  obj$data_info <- list("data_summary" = df,
                        "nca" = nca,
                        "nca_dose_norm" = nca_dose_norm,
                        "dose_norm_check" = dose_norm_check)

  obj$status <- status_data_info #data summarization complete

  return(obj)
}

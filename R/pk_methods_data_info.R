#' calculate summary data info
#'
#' Calculate summary data information, including non-compartmental analysis.
#'
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



  #get data summay
  if(obj$settings_preprocess$suppress.messages %in% FALSE){
    message("data_info.pk(): Getting data summary statistics\n")
  }
  data_summary_all <- data_summary(obj = obj,
                               newdata = NULL,
                               summary_group = vars(Chemical, Species))

  data_summary <- data_summary(obj = obj,
                               newdata = NULL,
                               summary_group = NULL)

  #perform NCA
  if(obj$settings_preprocess$suppress.messages %in% FALSE){
    message("data_info.pk(): Doing non-compartmental analysis\n")
  }
  nca_wide <- nca(obj = obj,
             newdata = NULL,
             nca_group = NULL,
             exclude = TRUE,
             dose_norm = FALSE)

  #do NCA dose-normalized
  # if(obj$settings_preprocess$suppress.messages %in% FALSE){
  #   message("data_info.pk(): Doing dose-normalized non-compartmental analysis\n")
  # }
  # nca_dose_norm <- nca(obj = obj,
  #                      newdata = NULL,
  #                      nca_group = vars(Chemical, Species, Route, Media),
  #                      exclude = TRUE,
  #                      dose_norm = TRUE)

  #get data flags:
  if(obj$settings_preprocess$suppress.messages %in% FALSE){
    message("data_info.pk(): Getting data flags\n")
  }

  #get grouping variables
  nca_group <- obj$settings_data_info$nca_group
  grp_vars <- sapply(nca_group,
                     rlang::as_label)


  df <- dplyr::inner_join(data_summary,
              nca_wide,
              by = grp_vars) %>%
    dplyr::mutate(`Cmax/Dose` = Cmax/Dose,
                  `AUC/Dose` = AUC_infinity/Dose) %>%
    dplyr::mutate(
    data_flag = ifelse(
      (n_obs - n_exclude) < 3,
      "Fewer than 3 non-excluded observations",
      NA_character_
    )
  ) %>% dplyr::mutate(
    data_flag = ifelse(
      n_detect < 3,
      paste2(data_flag,
             "Fewer than 3 non-excluded detected observations",
             sep = " | "),
      data_flag
    )
  ) %>% dplyr::mutate(
    data_flag = ifelse(
      Route %in% "oral" &
        (abs(tmax - tfirst_detect) < sqrt(.Machine$double.eps)) %in% TRUE &
                 !is.na(tmax),
      paste2(data_flag,
             "tmax is equal to time of first detect",
             sep = " | "),
      data_flag
    )) %>% dplyr::mutate(
      data_flag = ifelse(
        Route %in% "oral" &
          (abs(tmax - tlast_detect) < sqrt(.Machine$double.eps)) %in% TRUE &
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

  #Other data flags:
  #Check for obeying dose normalization by chemical, species, route, media

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

  obj$data_info <- list("data_summary_all" = data_summary_all,
    "data_summary_nca" = df,
                        "dose_norm_check" = dose_norm_check)

  obj$status <- status_data_info #data summarization complete

  return(obj)
}

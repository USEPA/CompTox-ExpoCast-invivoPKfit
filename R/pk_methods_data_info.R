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

  #get data summary
  if(obj$settings_preprocess$suppress.messages %in% FALSE){
    message("data_info.pk(): Getting data summary statistics\n")
  }

  #grouping for data summary:
  #data grouping, plus also Route and Media if not already included in data grouping
  summary_group <-   unique(
    c(obj$data_group,
      vars(Route, Media)
    )
  )

  data_summary_out <- data_summary(obj = obj,
                               newdata = NULL,
                               summary_group = summary_group
  )

  #do NCA dose-normalized
  if(obj$settings_preprocess$suppress.messages %in% FALSE){
    message("data_info.pk(): Doing dose-normalized non-compartmental analysis\n")
  }
  nca_dose_norm <- nca(obj = obj,
                       newdata = NULL,
                       nca_group = summary_group,
                       exclude = TRUE,
                       dose_norm = TRUE)
  #pivot wider
  #first get names of grouping vars
  grp_vars <- sapply(summary_group,
                     rlang::as_label)
#then pivot wider
  nca_dose_norm <-   nca_dose_norm %>%
    tidyr::pivot_wider(id_cols = tidyselect::all_of(grp_vars),
                       names_from = param_name,
                       values_from = param_value) %>%
    as.data.frame()

  #get data flags:
  if(obj$settings_preprocess$suppress.messages %in% FALSE){
    message("data_info.pk(): Getting data flags\n")
  }

  #get grouping variables

  grp_vars <- sapply(summary_group,
                     rlang::as_label)

  df <- dplyr::inner_join(data_summary_out,
              nca_dose_norm,
              by = grp_vars) %>%
   dplyr::mutate(
    data_flag = ifelse(
      Route %in% "oral" &
        (abs(tmax - tfirst_detect) < sqrt(.Machine$double.eps)) %in% TRUE &
                 !is.na(tmax),
    "tmax is equal to time of first detect",
      NA_real_
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
  #Check for obeying dose normalization by summary_group

  dose_norm_check <- do.call(dplyr::group_by,
                             args = c(list(df),
                             summary_group)) %>%
    dplyr::summarise(
      Cmax_fold_range = {
        tmprange <- suppressWarnings(range(`Cmax`, na.rm = TRUE))
        if(all(!is.finite(tmprange))) tmprange <- c(NA_real_, NA_real_)
        tmprange[2]/tmprange[1]
      },
      data_flag_Cmax = ifelse(Cmax_fold_range > 2,
                              "Cmax may not scale with dose. Cmax/Dose range > 2-fold across NCA groups for this Route/Media",
                              NA_character_),
      AUC_fold_range = {
        tmprange <- suppressWarnings(range(`AUC_infinity`, na.rm = TRUE))
        if(all(!is.finite(tmprange))) tmprange <- c(NA_real_, NA_real_)
        tmprange[2]/tmprange[1]
      },
      data_flag_AUC = ifelse(AUC_fold_range > 2,
                             "AUC_infinity may not scale with dose. AUC/Dose range > 2-fold across NCA groups for this Route/Media",
                             NA_character_),
    ) %>% as.data.frame()

  obj$data_info <- list("data_summary" = df,
                        "dose_norm_check" = dose_norm_check)

  obj$status <- status_data_info #data summarization complete

  return(obj)
}

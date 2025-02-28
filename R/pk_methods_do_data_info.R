#' calculate summary data info
#'
#' Calculate summary data information, including non-compartmental analysis.
#'
#'
#' @param obj An object of class [pk()].
#' @param ... Additional arguments. Not currently in use.
#' @return Object of class [pk()] with an added `$data_info` list containing
#' non-compartmental analysis results.
#' @export
#' @importFrom magrittr `%>%`
#' @author Caroline Ring
do_data_info.pk <- function(obj, ...) {

  # check status
  objname <- deparse(substitute(obj))
  status <- obj$status
  if (status >= status_data_info) {
    warning(objname,
            " current status is ",
            status,
            ". do_data_info.pk() will reset its status to ",
            status_data_info,
            ". Any results from later workflow stages will be lost.\n")
    # Here is where I need to implement logic for skipping if same data group
    prev_summary <- obj$data_info$data_summary %>%
      dplyr::select(Chemical:tlast_detect)
  }

  # if preprocessing not already done, do it
  if (obj$status < status_preprocess) {
    obj <- do_preprocess(obj)
  }

  if (obj$settings_preprocess$suppress.messages %in% FALSE) {
    message("do_data_info.pk(): Getting data summary info\n")
  }

  data <- obj$data

  # get data summary
  if (obj$settings_preprocess$suppress.messages %in% FALSE) {
    message("do_data_info.pk(): Getting data summary statistics\n")
  }

  # grouping for data summary:

  summary_group <- obj$settings_data_info$summary_group

  data_summary_out <- data_summary(obj = obj,
                                   newdata = NULL,
                                   summary_group = summary_group
  )
  # Here is step two of determining whether old data group is the same as "new" one
  if (status >= status_data_info) {
    id_summary <- identical(data_summary_out, prev_summary)
    if (id_summary) {
      message("Any changes made do not affect the outcome ",
              "of the non-compartmental analysis!\n")
      return(obj)
    }
  }

  # do NCA dose-normalized
  # for purposes of some data flags

  if (obj$settings_preprocess$suppress.messages %in% FALSE) {
    message("do_data_info.pk(): Doing dose-normalized non-compartmental analysis\n")
  }

  # if Dose, Route and Media not already in summary_group, add them for this
  grp_vars_summary <- sapply(summary_group,
                             rlang::as_label)

  if (all(c("Dose", "Route", "Media") %in% grp_vars_summary)) {
    nca_group <- summary_group

  } else {
    if (obj$settings_preprocess$suppress.messages %in% FALSE) {
      message(paste0("Grouping for NCA must include Dose, Route, and Media, but ",
                     "summary_group is: ",
                     toString(grp_vars_summary),
                     "\n",
                     "Dose, Route, and/or Media are being added to the grouping for purposes of NCA."))
    }
    nca_group <- union(summary_group,
                       vars(Dose, Route, Media))
  }

  nca_dose_norm_long <- nca(obj = obj,
                       newdata = NULL,
                       nca_group = nca_group,
                       exclude = TRUE,
                       dose_norm = TRUE)

  # pivot wider
  grp_vars_nca <- sapply(nca_group,
                         rlang::as_label)

  nca_dose_norm <- nca_dose_norm_long %>%
    tidyr::pivot_wider(id_cols = tidyselect::all_of(c(grp_vars_nca,
                                                      "dose_norm")),
                       names_from = param_name,
                       values_from = param_value) %>%
    as.data.frame()

  # get data flags:
  if (obj$settings_preprocess$suppress.messages %in% FALSE) {
    message("do_data_info.pk(): Getting data flags based on NCA\n")
  }

  # get summary for nca_group if it is different from summary_group
if (length(setdiff(nca_group, summary_group)) > 0 |
   length(setdiff(summary_group, nca_group)) > 0) {
  data_summary_nca <- data_summary(obj = obj,
                                   newdata = NULL,
                                   summary_group = nca_group
  )
} else {
  data_summary_nca <- data_summary_out
}


  # Assess various flags
  df <- dplyr::inner_join(data_summary_nca,
                          nca_dose_norm,
                          by = grp_vars_nca) %>%
    dplyr::mutate(
      # Is tmax for oral data the first or last detected timepoint?
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
          # CLtot/Fgutabs must be positive and AUC_infinity must also be positive
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
        )

  # Other data flags:
  # Check for obeying dose normalization by summary_group - Dose
  dose_norm_check <- df %>%
    dplyr::group_by(!!!summary_group) %>%
    dplyr::ungroup(Dose) %>%
    dplyr::reframe(
      n_dose_groups = dplyr::n_distinct(Dose),
      Cmax_fold_range = {
        tmprange <- suppressWarnings(range(`Cmax`, na.rm = TRUE))
        if (all(!is.finite(tmprange))) tmprange <- c(NA_real_, NA_real_)
        tmprange[2] / tmprange[1]
      },
      data_flag_Cmax = ifelse(Cmax_fold_range > 2,
                              "Cmax may not scale with dose. Cmax/Dose range > 2-fold across dose groups for this group",
                              NA_character_),
      AUC_fold_range = {
        tmprange <- suppressWarnings(range(`AUC_infinity`, na.rm = TRUE))
        if (all(!is.finite(tmprange))) tmprange <- c(NA_real_, NA_real_)
        tmprange[2] / tmprange[1]
      },
      data_flag_AUC = ifelse(AUC_fold_range > 2,
                             "AUC_infinity may not scale with dose. AUC/Dose range > 2-fold across dose groups for this group",
                             NA_character_)) %>%
    as.data.frame()

  obj$data_info <- list("data_summary" = data_summary_out,
                        "data_flags" = df,
                        "dose_norm_check" = dose_norm_check,
                        "nca" = nca_dose_norm_long)

  obj$status <- status_data_info # data summarization complete

  return(obj)
}

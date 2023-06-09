#' NCA for a `pk` object
#'
#' Non-compartmental analysis for a `pk` object
#'
#' Perform non-compartmental analysis of data in a `pk` object (or optionally,
#' new data), using data groupings defined by `get_nca_group()` for the `pk`
#' object (or optionally, new groupings). If you provide both `newdata` and
#' `nca_group`, then everything in the `pk` object will be ignored and you will
#' simply be doing NCA *de novo* (which may be what you want).
#'
#' @param obj A [pk()] model object. Must be fitted, or the function will exit
#'   with an error.
#' @param newdata Optional: A `data.frame` containing new data for which to
#'   compute the TK stats. Must contain at least variables `Chemical`,
#'   `Species`, `Route`, `Dose`, `Conc`, `Dose.Units`, `Conc.Units`, either
#'   `Time_trans.Units` or `Time.Units`, and any other variables named in
#'   `tk_grouping`. Default `NULL`, to use the data in `get_data(obj)`.
#' @param nca_group A list of variables provided using a `dplyr::vars()` call.
#'   The data (either `newdata` or `obj$data`) will be grouped according to the
#'   unique combinations of these variables. For each unique combination of
#'   these variables in the data, a set of TK statistics will be computed. The
#'   default is `NULL`, to use the same data grouping that was set in
#'   [stat_nca()] for the `pk` object. However, you may specify a different data
#'   grouping if you wish.
#' @param exclude Logical: `TRUE` to group the data for NCA after removing any
#'   observations in the data marked for exclusion (if there is a variable
#'   `exclude` in the data, an observation is marked for exclusion when `exclude
#'   %in% TRUE`). `FALSE` to include all observations when grouping the data for
#'   NCA, regardless of exclusion status. Default `TRUE`.
#' @param dose_norm Logical: `TRUE` to perform NCA after dose-normalizing
#'   concentrations. `FALSE` (default) to perform NCA on un-transformed
#'   concentrations.
#' @return A `data.frame` with variables including all the grouping variables in
#'   `nca_group`, `nca_group_id`; `design` (the auto-detected study design for
#'   this group); `param_name` (the name of the NCA parameter); `param_value`
#'   (the NCA parameter value); `param_sd_z` (standard deviation of the
#'   estimated NCA parameter value, if avaialble); `param_units` (the units of
#'   the NCA parameter, derived from the units of the data).
#' @export
#' @author Caroline Ring

nca.pk <- function(obj,
                   newdata = NULL,
                   nca_group = NULL,
                   exclude = TRUE,
                   dose_norm = FALSE){

  suppress.messages <- obj$settings_preprocess$suppress.messages

  if(is.null(nca_group)){
    nca_group <- obj$settings_data_info$nca_group
  }

  if(is.null(newdata)) newdata <- obj$data

  grp_vars <- sapply(nca_group,
                     rlang::as_label)

  newdata_ok <- check_newdata(newdata = newdata,
                              olddata = obj$data,
                              req_vars = union(
                                c("Chemical",
                                  "Species",
                                  "Time",
                                  "Time.Units",
                                  "Dose",
                                  "Conc",
                                  "Dose.Units",
                                  "Conc.Units",
                                  "Route"),
                                grp_vars),
                              exclude = exclude)

  #if exclude = FALSE, then treat all observations as included
  if(exclude %in% FALSE){
    newdata$exclude <- FALSE
  }

  if(dose_norm %in% TRUE){
    newdata$Conc <- newdata$Conc/newdata$Dose
    newdata$Dose <- newdata$Dose/newdata$Dose
    newdata$Conc.units <- paste0(newdata$Conc.Units,
                                 "/",
                                 newdata$Dose.Units)
  }

  #do NCA

  nca_out <- do.call(dplyr::group_by,
          args =c(list(newdata),
                  nca_group)) %>%
    dplyr::summarise(Conc.Units = unique(Conc.Units),
                     Time_trans.Units = unique(Time_trans.Units),
                     Dose.Units = unique(Dose.Units),
                     {
                       if(suppress.messages %in% FALSE){
                         cur_data_summary <- dplyr::inner_join(get_data_summary(obj,
                                                                                summary_group = nca_group),
                                                               dplyr::cur_group(),
                                                               by = grp_vars) %>%
                           as.data.frame
                         message(paste("nca.pk(): Doing",
                                       ifelse(dose_norm %in% TRUE,
                                              "dose-normalized",
                                              "non-dose-normalized"),
                                       "NCA for the following data:"))
                         print(cur_data_summary)
                       }
                       calc_nca(time = Time_trans[exclude %in% FALSE], #calculate NCA
                              dose = Dose[exclude %in% FALSE],
                              conc = Conc[exclude %in% FALSE],
                              detect = Detect[exclude %in% FALSE],
                              route = unique(Route[exclude %in% FALSE]),
                              series_id = Series_ID[exclude %in% FALSE])
                       }
                     ) %>%
    dplyr::mutate(param_units = dplyr::case_when( #derive NCA param units from data units
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
                     Dose.Units))



  return(nca_out)
}

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
#'   `Species`, `Route`, `Dose`, `Conc`, `Dose.Units`, `Conc.Units`, and
#'   `Time.Units`, and any other variables named in
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
#' @param suppress.messages Logical.
#' @param ... Additional arguments. Currently not in use.
#' @return A `data.frame` with variables including all the grouping variables in
#'   `nca_group`, `nca_group_id`; `design` (the auto-detected study design for
#'   this group); `param_name` (the name of the NCA parameter); `param_value`
#'   (the NCA parameter value); `param_sd_z` (standard deviation of the
#'   estimated NCA parameter value, if available); `param_units` (the units of
#'   the NCA parameter, derived from the units of the data).
#' @export
#' @author Caroline Ring

nca.pk <- function(obj,
                   newdata = NULL,
                   nca_group = NULL,
                   exclude = TRUE,
                   dose_norm = FALSE,
                   suppress.messages = NULL,
                   ...) {

  if (is.null(suppress.messages)) {
    suppress.messages <- obj$settings_preprocess$suppress.messages
  }

  if (is.null(nca_group)) {
    nca_group <- obj$settings_data_info$summary_group
  }

  if (is.null(newdata)) newdata <- obj$data

  grp_vars <- sapply(nca_group,
                     rlang::as_label)

  # Create a new dose variable to handle dose-normalization or not
  # If dose-normalized, grouping need not include original Dose column

  if (dose_norm %in% TRUE) {
    # check_nca_group to ensure it includes  Route, and Media
    # since we can only do NCA for a single Route, and a single Media at a time
    # (when dose_norm is TRUE, Dose = 1, so we already have a single Dose)
    if (!(all(c("Route",
               "Media") %in%
             grp_vars))) {
      stop(paste0("nca.pk(): When dose_norm == TRUE, nca_group must include all of Route, Media.\n",
                  "nca_group is: ",
                  toString(grp_vars)
      ))
    }
  } else {
    # check_nca_group to ensure it includes Dose, Route, and Media
    # since we can only do NCA for a single Dose, a single Route, and a single Media at a time
    if (!(all(c("Dose",
               "Route",
               "Media") %in%
             grp_vars))) {
      stop(paste0("When dose_norm == FALSE, nca_group must include all of Dose, Route, Media.\n",
                 "nca.pk(): nca_group is: ",
                 toString(grp_vars)
      ))
    }
  }

  if (dose_norm %in% TRUE) {
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
  } else { # if dose_norm == TRUE
    # ignore Dose column in check because it is technically not required
    newdata_ok <- check_newdata(newdata = newdata,
                                olddata = obj$data,
                                req_vars = union(
                                  c("Chemical",
                                    "Species",
                                    "Time",
                                    "Time.Units",
                                    "Conc",
                                    "Dose.Units",
                                    "Conc.Units",
                                    "Route"),
                                  grp_vars),
                                exclude = exclude)
  }


  # if exclude = FALSE, then treat all observations as included
  if (exclude %in% FALSE) {
    newdata$exclude <- FALSE
  }

  if (dose_norm %in% TRUE) {
    newdata$Conc_nca <- newdata$Conc / newdata$Dose
    newdata$Dose_nca <- 1.0
    newdata$Conc_nca.Units <- paste0("(",
                                     newdata$Conc.Units,
                                     ")",
                                 "/",
                                 "(",
                                 newdata$Dose.Units,
                                 ")")
  } else {
    newdata$Conc_nca <- newdata$Conc
    newdata$Dose_nca <- newdata$Dose
    newdata$Conc_nca.Units <- newdata$Conc.Units
  }

  if (suppress.messages %in% FALSE) {
    message(paste("nca.pk(): Doing",
                                  ifelse(dose_norm %in% TRUE,
                                         "dose-normalized",
                                         "non-dose-normalized"),
                                  "NCA by the following grouping:",
                  toString(grp_vars)
                  )
            )
  }

  # do NCA
    nca_out <- newdata %>%
      dplyr::group_by(!!!nca_group) %>%
      dplyr::reframe(Conc.Units = unique(Conc_nca.Units),
                     Time.Units = unique(Time.Units),
                     Dose.Units = unique(Dose.Units),
                     dose_norm = dose_norm,
                     {
                       calc_nca(time = Time[exclude %in% FALSE], # calculate NCA
                                dose = Dose_nca[exclude %in% FALSE],
                                conc = Conc_nca[exclude %in% FALSE],
                                detect = Detect[exclude %in% FALSE],
                                route = unique(Route[exclude %in% FALSE]),
                                series_id = Series_ID[exclude %in% FALSE])
                     }
      ) %>%
      dplyr::mutate(
        param_units = dplyr::case_when( # derive NCA param units from data units
          param_name %in% c("AUC_tlast",
                            "AUC_infinity") ~ paste0("(",
                                                    Conc.Units,
                                                    ")",
                                                    "*",
                                                    Time.Units),
          param_name %in% "AUMC_infinity" ~ paste0("(",
                                                  Conc.Units,
                                                  ")",
                                                  "*",
                                                  Time.Units,
                                                  "*",
                                                  Time.Units),
          param_name %in% c("MRT",
                            "MTT",
                            "halflife",
                            "tmax") ~ Time.Units,
          param_name %in% c("CLtot",
                            "CLtot/Fgutabs") ~ paste0("L/",
                                                      Time.Units),
          param_name %in% "Vss" ~ paste0("(",
                                         Conc.Units,
                                         ")",
                                         "/",
                                         "(",
                                         Dose.Units,
                                         ")"),
          param_name %in% "Cmax" ~ Conc.Units
        )) %>%
      dplyr::select(-c(Conc.Units,
                       Time.Units,
                       Dose.Units))

  return(nca_out)
}

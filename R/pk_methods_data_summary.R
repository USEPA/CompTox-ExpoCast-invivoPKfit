#' Data summary for a `pk` object
#'
#' Calculate data summary statistics for a `pk` object
#'
#' Get summary statistics for data in a `pk` object (or optionally, new data),
#' using data groupings defined by `get_nca_group()` for the `pk` object (or
#' optionally, new groupings). If you provide both `newdata` and
#' `summary_group`, then everything in the `pk` object will be ignored and you
#' will simply be doing data summary *de novo* (which may be what you want).
#'
#' Summary statistics include, for each group:
#' \itemize{
#' \item `n_obs`: the number of observations
#' \item `n_exclude`: The number of excluded observations
#' \item `n_detect`: The number of non-excluded detected observations
#' \item `n_series_id`: The number of unique series IDs
#' \item `n_timepts`: The number of unique time points
#' \item `n_ref`: The number of unique reference IDs
#' \item `tlast`: The time of the latest non-excluded observation
#' \item `tlast_detect`: The time of the latest non-excluded detected observation
#' \item `tfirst`: The time of the earliest non-excluded observation
#' \item `tfirst_detect`: The time of the earliest non-excluded detected observation
#' }
#'
#' @param obj A [pk()] model object. Must be fitted, or the function will exit
#'   with an error.
#' @param newdata Optional: A `data.frame` containing new data for which to
#'   compute the TK stats. Must contain at least variables `Chemical`,
#'   `Species`, `Route`, `Dose`, `Conc`, `Dose.Units`, `Conc.Units`, either
#'   `Time_trans.Units` or `Time.Units`, and any other variables named in
#'   `tk_grouping`. Default `NULL`, to use the data in `get_data(obj)`.
#' @param summary_group A list of variables provided using a `dplyr::vars()`
#'   call. The data (either `newdata` or `obj$data`) will be grouped according
#'   to the unique combinations of these variables. For each unique combination
#'   of these variables in the data, a set of summary statistics will be
#'   computed. The default is `NULL`, to use the same data grouping that was set
#'   in [stat_nca_group()] for the `pk` object. However, you may specify a different
#'   data grouping if you wish.
#' @param ... Additional arguments. Not in use.
#' @return A `data.frame` with variables including all the grouping variables in
#'   `summary_group`, `group_id`; `param_name` (the name of the summary
#'   statistic; see Details); `param_value` (the summary statistic value);  `param_units`
#'   (the units of the summary statistic, derived from the units of the data).
#' @export
#' @author Caroline Ring, Gilberto Padilla Mercado

data_summary.pk <- function(obj,
                            newdata = NULL,
                            summary_group = NULL,
                            ...) {

  if (is.null(summary_group)) {
    summary_group <- get_nca_group.pk(obj)
  }

  if (is.null(newdata)) newdata <- get_data.pk(obj)

  grp_vars <- sapply(summary_group,
                     rlang::as_label)

  newdata_ok <- check_newdata(newdata = newdata,
                              olddata = obj$data,
                              req_vars = union(
                                c("Chemical",
                                  "Species",
                                  "Time",
                                  "Time.Units",
                                  "Dose",
                                  "Dose.Units",
                                  "Conc.Units",
                                  "Route",
                                  "Media"),
                                grp_vars),
                              exclude = TRUE)


  data_summary <- dplyr::group_by(newdata, !!!summary_group) |>
    dplyr::reframe(
      n_obs = dplyr::n(), # summary stats on data
      n_exclude = sum(exclude %in% TRUE),
      n_detect = sum(Detect %in% TRUE & exclude %in% FALSE),
      n_series_id = length(unique(Series_ID[exclude %in% FALSE])),
      n_study_id = length(unique(Study_ID[exclude %in% FALSE])),
      n_timepts = length(unique(Time_trans[exclude %in% FALSE])),
      n_ref = length(unique(Reference[exclude %in% FALSE])),
      tfirst = ifelse(any(exclude %in% FALSE),
                      min(Time_trans[exclude %in% FALSE]),
                      NA_real_),
      tfirst_detect = ifelse(any(Detect %in% TRUE & exclude %in% FALSE),
                             min(Time_trans[Detect %in% TRUE & exclude %in% FALSE]),
                             NA_real_),
      tlast = ifelse(any(exclude %in% FALSE),
                     max(Time_trans[exclude %in% FALSE]),
                     NA_real_),
      tlast_detect = ifelse(any(Detect %in% TRUE & exclude %in% FALSE),
                            max(Time_trans[Detect %in% TRUE & exclude %in% FALSE]),
                            NA_real_)

    )

  data_summary <- as.data.frame(data_summary)

  return(data_summary)
}

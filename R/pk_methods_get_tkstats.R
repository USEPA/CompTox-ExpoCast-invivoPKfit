#' Get TK stats
#'
#' Extract derived TK statistics from a fitted [pk()] model object.
#'
#' After fitting model parameters (e.g. elimination rate, volume of
#' distribution, absorption rate, bioavailability), it can be useful to derive
#' summary toxicokinetic statistics such as total clearance rate, half-life,
#' peak concentration, AUC_inf (the area under the concentration-time curve when
#' time goes to infinity), etc.
#'
#' Many of these TK statistics depend not only on chemical and species, but also
#' on route, media (tissue), and dose. Therefore, TK stats need to be computed
#' for a specific set of Chemical, Species, Route, Media, and Dose.
#'
#' TK statistics for a defined [pk_model()] object are computed using the
#' function named in the model's `tkstats_fun`. For the built-in models, the
#' `tkstats_fun` functions are the following. See the documentation for the
#' individual functions for details on what TK stats are calculated for each
#' model, and how they are calculated.

#' -`model_1comp`: [tkstats_1comp()]
#' -`model_2comp`: [tkstats_2comp()]
#' -`model_flat`: [tkstats_flat()]
#'
#' @param obj A [pk()] model object. Must be fitted, or the function will exit
#'   with an error.
#' @param newdata Optional: A `data.frame` containing new data for which to
#'   compute the TK stats. Must contain at least variables `Chemical`,
#'   `Species`, `Route`, `Media`, `Dose`, `Dose.Units`, `Conc.Units`, either
#'   `Time_trans.Units` or `Time.Units`, and any other variables named in
#'   `tk_grouping`. Default `NULL`, to use the data in `obj$data`.
#' @param tk_group A list of variables provided using a `dplyr::vars()` call.
#'   The data (either `newdata` or `obj$data`) will be grouped according to the
#'   unique combinations of these variables. For each unique combination of
#'   these variables in the data, a set of TK statistics will be computed. The
#'   default is `obj$settings_data_info$nca_group`, to derive TK statistics for
#'   the same groups of data as non-compartmental analysis statistics. With the
#'   default, you can directly compare e.g. a model-predicted AUC_inf to the
#'   corresponding NCA-estimated AUC_inf. However, you may specify a different
#'   data grouping if you wish. Each group should have a unique combination of
#'   `Chemical`, `Species`, `Route`, `Media`, and `Dose`, because the TK stats
#'   depend on these values, and it is required to have one unique set of TK
#'   stats per group.
#' @param model Character: One or more of the models fitted. Default `NULL` to
#'   return TK stats for all models.
#' @param method Character: One or more of the [optimx::optimx()] methods used.
#'   Default `NULL` to return TK stats for all methods.
#' @param exclude Logical: `TRUE` to get the TK groupings after removing any
#'   observations in the data marked for exclusion (if there is a variable
#'   `exclude` in the data, an observation is marked for exclusion when  `TRUE`).
#'   `FALSE` to include all observations when getting the TK
#'   groupings, regardless of exclusion status. Default `TRUE`.
#' @param vol_unit Character: Specifies the unit of volume. Defaults to "L" for liters.
#' @param dose_norm Logical: `TRUE` (default) specifies whether the concentrations are dose-normalized.
#' @param ... Additional arguments not currently in use.
#' @return  A data.frame with one row for each `data_group`, `model` and `method`
#'   with the variables in the `data.frame` returned by the `tkstats_fun` for
#'   its corresponding model.
#'   (For the built-in models `model_flat`, `model_1comp`, and `model_2comp`, these
#'   variables are `param_name` and `param_value`.)
#' @export
#' @author Caroline Ring, Gilberto Padilla Mercado
#' @family methods for fitted pk objects
get_tkstats.pk <- function(obj,
                           newdata = NULL,
                           tk_group = NULL,
                           model = NULL,
                           method = NULL,
                           exclude = TRUE,
                           vol_unit = "L",
                           dose_norm = TRUE, ...){

  #ensure that the model has been fitted
  check <- check_required_status(obj = obj,
                                 required_status = status_fit)
  if (!(check %in% TRUE)) {
    stop(attr(check, "msg"))
  }

  if (is.null(model)) model <- names(obj$stat_model)
  if (is.null(method)) method <- obj$settings_optimx$method
  if (is.null(newdata)) newdata <- obj$data
  if (is.null(tk_group)) tk_group <- obj$settings_data_info$nca_group

  method_ok <- check_method(obj = obj, method = method)
  model_ok <- check_model(obj = obj, model = model)

  grp_vars <- sapply(tk_group,
                     rlang::as_label)
  data_grp_vars <- sapply(obj$data_group, rlang::as_label)

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
                              exclude = exclude)

  #if exclude = TRUE, remove excluded observations
  if (exclude %in% TRUE) {
    newdata <- subset(newdata, exclude %in% FALSE)
  }


  #check that tk_group is valid: it must produce groups with a unique
  #combination of Chemical, Species, Route, Media, and Dose
  newdata_grouped <- newdata %>%
    dplyr::group_by(!!!tk_group) %>%
    dplyr::distinct(Chemical, #for each tk_group, take distinct rows by these variables
                    Species,
                    Route,
                    Media,
                    Dose,
                    Time.Units,
                    Dose.Units,
                    Conc.Units)

  newdata_grouped_count <- newdata_grouped %>%
    dplyr::count(name = "N") %>%
    dplyr::ungroup() #how many distinct rows per group?


  #if more than one distinct row per group, stop
  if (any(newdata_grouped_count$N > 1)) {
    stop("tk_group does not produce groups with unique combinations of Chemical, Species, Route, Media, and Dose.")
  }

  all_coefs <- coef(obj,
                    model = model,
                    method = method,
                    drop_sigma = TRUE)

  model_df <- data.frame(model = sapply(obj$stat_model, `[[`, "name"),
                         tk_fun = sapply(obj$stat_model, `[[`, "tkstats_fun"))

  tkstats_all <- dplyr::left_join(all_coefs,
                                  newdata_grouped,
                                  by = union(data_grp_vars, "Time.Units"),
                                  relationship = "many-to-many") %>%
    dplyr::left_join(model_df,
                     by = "model",
                     relationship = "many-to-many") %>%
    dplyr::group_by(!!!obj$data_group, model, method,
                    coefs_vector, tk_fun) %>%
    tidyr::nest(.key = "c_data")

  tkstats_all <- tkstats_all %>%
    dplyr::reframe(tkstats_df = purrr::map(c_data,
                                    \(x) {
                                      x %>%
                                        dplyr::rowwise() %>%
                                        dplyr::mutate(
                                          TKstats = sapply(
                                            coefs_vector,
                                            FUN = tk_fun,
                                            route = Route,
                                            medium = Media,
                                            dose = ifelse(dose_norm == TRUE, 1, Dose),
                                            time_unit = "hours", # coef() is standardized now
                                            conc_unit = Conc.Units,
                                            vol_unit = "L",
                                            simplify = FALSE,
                                            USE.NAMES = TRUE
                                          )
                                        ) %>% dplyr::ungroup()
                                    })) %>%
    tidyr::unnest(cols = c(tkstats_df))

  tkstats_all <- tkstats_all %>%
    dplyr::mutate(TKstats_wide  = purrr::map(TKstats,
                               \(x) {
                                 dplyr::select(x, -param_units) %>%
                                   tidyr::pivot_wider(names_from = param_name,
                                                      values_from = param_value)
                               })) %>%
    tidyr::unnest(cols = c(TKstats_wide)) %>%
    dplyr::group_by(!!!tk_group)


  # Final filtering of tkstats_all

  tkstats_all <- tkstats_all[!names(tkstats_all) %in% c("coefs_vector", "tk_fun", "TKstats")] %>%
    dplyr::relocate(!!!tk_group, Dose.Units, Conc.Units) %>%
    dplyr::ungroup()

  if (dose_norm) {
    tkstats_all <- tkstats_all %>%
      dplyr::select(!Dose)

    message("Dose column removed because these TK statistics are dose normalized")
  }

  tkstats_all <- tkstats_all %>% dplyr::distinct()

  return(tkstats_all)
}

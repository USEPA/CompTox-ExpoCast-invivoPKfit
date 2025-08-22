#' Data preprocessing settings
#'
#' @param routes_keep Character: A list of routes to keep. Data will be filtered
#'   so that the harmonized variable `Route` includes only values in
#'   `routes_keep`. Default is `c("oral", "iv")`.
#' @param media_keep Character: A list of media to keep. Data will be filtered
#'   so that the harmonized variable name `Media` includes only values in
#'   `media_keep`. Default is `c("blood", "plasma")`.
#' @param ratio_conc_dose Numeric: The ratio of mass units of observed
#'   concentrations to mass units of applied doses. Default 1, to indicate the
#'   same mass units are used for both.
#' @param impute_loq TRUE or FALSE: Whether to impute missing LOQ values.
#' @param calc_loq_factor  A numeric factor used for imputing missing LOQ.
#'   Within each group defined in `loq_group`, any missing LOQ values will be
#'   imputed as the minimum detected Value in the group, multiplied by
#'   `calc_loq_factor`. Default 0.9.
#' @param impute_sd TRUE or FALSE: Whether to impute missing SD values.
#' @param keep_data_original TRUE or FALSE: Whether to keep original data after pre-processing.
#' @param suppress.messages TRUE or FALSE: Whether to suppress verbose messages.
#'   Default FALSE.
#' @param ... Any additional arguments. Currently ignored.
#'
#' @return An object of class `pk_settings_preprocess`. This is a named list of
#'   the arguments provided to this function and their values.
#' @author Caroline Ring
#' @export
settings_preprocess <- function(routes_keep = c("oral", "iv"),
                                media_keep = c("blood", "plasma", "excreta"),
                                ratio_conc_dose = 1,
                                impute_loq = TRUE,
                                calc_loq_factor = 0.9,
                                impute_sd = TRUE,
                                keep_data_original = TRUE,
                                suppress.messages = FALSE,
                                ...) {
  # get arguments and values
  argg <- c(as.list(environment()), list(...))
  this_settings_preprocess <- argg
  # set class
  class(this_settings_preprocess) <- c(class(this_settings_preprocess), "pkproto", "pk_settings_preprocess")

  return(this_settings_preprocess)
}

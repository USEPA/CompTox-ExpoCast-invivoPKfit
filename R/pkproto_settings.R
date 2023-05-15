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
#' @param loq_group A list of variables, specified using a call to
#'   [dplyr::vars()]. These should be harmonized variable names. Unique
#'   combinations of these variables define groups of data. Within each group,
#'   any missing LOQ values will be imputed as the minimum detected Value in the
#'   group, multiplied by `calc_loq_factor`. Default is `dplyr::vars(Chemical,
#'   Species, Reference, Media)`.
#' @param calc_loq_factor  A numeric factor used for imputing missing LOQ.
#'   Within each group defined in `loq_group`, any missing LOQ values will be
#'   imputed as the minimum detected Value in the group, multiplied by
#'   `calc_loq_factor`. Default 0.45.
#' @param impute_sd TRUE or FALSE: Whether to impute missing SD values.
#' @param sd_group A list of variables, specified using a call to
#'   [dplyr::vars()]. These should be harmonized variable names. Unique
#'   combinations of these variables define groups of data. Within each group,
#'   any missing SD values will be imputed as the minimum non-missing SD value
#'   in the group. Default is `dplyr::vars(Chemical, Species, Reference,
#'   Media)`.
#' @param suppress.messages TRUE or FALSE: Whether to suppress verbose messages.
#'   Default FALSE.
#' @param ... Any additional arguments. Currently ignored.
#'
#' @return An object of class `pk_settings_preprocess`. This is a named list of the
#'   arguments provided to this function and their values.
#' @author Caroline Ring
#' @export
settings_preprocess <- function(routes_keep = c("oral", "iv"),
                             media_keep = c("blood", "plasma"),
                          ratio_conc_dose = 1,
                             impute_loq = TRUE,
                          loq_group = dplyr::vars(Chemical, Species, Reference, Media),
                          calc_loq_factor = 0.45,
                             impute_sd = TRUE,
                          sd_group = dplyr::vars(Chemical, Species, Reference, Media),

                             suppress.messages = FALSE,
                          ...){
#get arguments and values
argg <- c(as.list(environment()), list(...))
this_settings_preprocess <- argg
#set class
class(this_settings_preprocess) <- c(class(this_settings_preprocess), "pkproto", "pk_settings_preprocess")

return(this_settings_preprocess)
}

#' Data info settings
#'
#' @param nca_group A set of variables. Data will be split into groups according
#'   to unique combinations of these variables, and non-compartmental analysis
#'   will be performed separately on each group. Default `dplyr::vars(Chemical,
#'   Species, Reference, Route, Media, Dose)`.
#' @param ... Other arguments (currently ignored).
#' @return An object of class `c(pkproto, pk_settings_data_info)`
#' @export
settings_data_info <- function(nca_group = dplyr::vars(Chemical, Species, Reference, Route, Media, Dose),
                               ...){
  #get arguments and values
  argg <- c(as.list(environment()), list(...))
  this_settings_data_info <- argg
  #set class
  class(this_settings_data_info) <- c(class(this_settings_data_info), "pkproto", "pk_settings_data_info")

  return(this_settings_data_info)
}

#' `optimx` optimizer settings
#'
#' @param method The name(s) of optimization methods to be used. See
#'   [optimx::optimx()] for options. Default is `"bobyqa"`.
#' @param itnmax The maximum number of iterations; as in [optimx::optimx()].
#' @param hessian Whether to compute the Hessian at the final set of parameters; as in [optimx::optimx()].
#' @param control A list of control parameters for the optimizer; see [optimx::optimx()] for options and details.
#' @return An object of class `pk_settings`.
#' @author Caroline Ring
settings_optimx <- function(method = c("bobyqa", "L-BFGS-B"),
                               itnmax = 1e6,
                               hessian = FALSE,
                               control = list(kkt = FALSE),
                            ...){
  #get arguments and values
  argg <- c(as.list(environment()), list(...))
  this_settings_optimx <- argg
  #set class
  class(this_settings_optimx) <- c(class(this_settings_optimx), "pkproto", "pk_settings_optimx")

  return(this_settings_optimx)
}

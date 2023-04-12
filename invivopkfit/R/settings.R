#' Data preprocessing settings
#'
#' @param ratio_conc_to_dose
#' @param calc_loq_factor
#' @param routes_keep
#' @param media_keep
#' @param impute_loq
#' @param impute_sd
#' @param suppress.messages
#'
#' @return An object of class `pk_settings`.
#' @author Caroline Ring
settings_data.pk <- function(ratio_conc_to_dose = 1,
                             calc_loq_factor = 0.45,
                             routes_keep = c("po", "iv"),
                             media_keep = c("blood", "plasma"),
                             impute_loq = TRUE,
                             impute_sd = TRUE,
                             suppress.messages = FALSE){
settings_data <- list(name = "data_settings")
#get arguments and values
argg <- c(as.list(environment()), list(...))
settings_data$value <- argg
#set class
class(settings_data) <- c(class(settings_data), "pkproto", "pk_settings")

return(settings_data)
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
settings_optimx.pk <- function(method = "bobyqa",
                               itnmax = 1e6,
                               hessian = FALSE,
                               control = list(kkt = FALSE)){
  settings_optimx <- list(name = "optimx_settings")
  #get arguments and values
  argg <- c(as.list(environment()), list(...))
  settings_optimx$value <- argg
  #set class
  class(settings_optimx) <- c(class(settings_optimx), "pkproto", "pk_settings")

  return(settings_optimx)
}

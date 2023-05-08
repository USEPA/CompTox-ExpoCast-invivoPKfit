#' Data preprocessing settings
#'

#' @param routes_keep
#' @param media_keep
#' @param impute_loq
#' @param loq_group
#' @param calc_loq_factor
#' @param impute_sd
#' @param suppress.messages
#'
#' @return An object of class `pk_data_settings`.
#' @author Caroline Ring
settings_data <- function(routes_keep = c("oral", "iv"),
                             media_keep = c("blood", "plasma"),
                             impute_loq = TRUE,
                          loq_group = vars(Chemical, Species, Reference, Media),
                          calc_loq_factor = 0.45,
                             impute_sd = TRUE,
                          sd_group = vars(Chemical, Species, Reference, Media),
                             suppress.messages = FALSE,
                          ...){
#get arguments and values
argg <- c(as.list(environment()), list(...))
this_settings_data <- argg
#set class
class(this_settings_data) <- c(class(this_settings_data), "pkproto", "pk_data_settings")

return(this_settings_data)
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
settings_optimx <- function(method = "bobyqa",
                               itnmax = 1e6,
                               hessian = FALSE,
                               control = list(kkt = FALSE),
                            ...){
  #get arguments and values
  argg <- c(as.list(environment()), list(...))
  this_settings_optimx <- argg
  #set class
  class(this_settings_optimx) <- c(class(this_settings_optimx), "pkproto", "pk_optimx_settings")

  return(this_settings_optimx)
}

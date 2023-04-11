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
settings_data <- list(name = "data")
#get arguments and values
argg <- c(as.list(environment()), list(...))
settings_data$value <- argg
#set class
class(settings_data) <- c(class(settings_data), "pkproto", "pk_settings")

return(settings_data)
}

#' `optimx` optimizer settings
#'
#'
settings_optimx.pk <- function(method = "bobyqa",
                               itnmax = 1e6,
                               hessian = FALSE,
                               control = list(kkt = FALSE)){
  settings_optimx <- list(name = "optimx")
  #get arguments and values
  argg <- c(as.list(environment()), list(...))
  settings_optimx$value <- argg
  #set class
  class(settings_optimx) <- c(class(settings_optimx), "pkproto", "pk_settings")

  return(settings_optimx)
}

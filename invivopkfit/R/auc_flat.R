#' AUC for flat model
#'
#'@param time A vector of times in hours
#'@param params A named list of model parameter values. Must include:
#'\describe{
#'\item{A}{Average concentration}
#'}
#'@param dose Dose in mg/kg
#'@param iv.dose TRUE for single IV bolus dose; FALSE for single oral dose
#'
#'@return A vector of AUC values evaluated at each time point in \code{time}.
#'
#'@author Caroline Ring, John Wambaugh, Chris Cook
#'
#' @export cp_flat
auc_flat <- function(params, time, dose, iv.dose){
  auc <- params$A * time
  auc[auc<0] <- 0
  auc[!is.finite(auc)] <- 0
  return(auc)
}

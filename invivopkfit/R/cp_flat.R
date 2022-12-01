#'Flat model
#'
#'@param params A named list of model parameter values. Must include: \describe{
#'  \item{A}{Average concentration} }
#'@param time A numeric vector of times in hours.
#'@param dose A numeric vector of doses in mg/kg
#'@param iv.dose A logical vector: TRUE for single IV bolus dose; FALSE for
#'  single oral dose. Not used, but must be present for compatibility with other
#'  model functions.
#'
#'@return A vector of plasma concentration values (mg/L) corresponding to \code{time}.
#'
#'@author Caroline Ring, John Wambaugh, Chris Cook
#'
#'@export cp_flat

cp_flat <- function(time, params, dose, iv.dose) {

  cp <- params$A * dose

  return(cp)
}

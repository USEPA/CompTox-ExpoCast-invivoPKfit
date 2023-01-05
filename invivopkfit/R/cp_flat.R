#'Flat model
#'
#'Evaluates a "flat" model for concentration vs. time: `Cp = A * dose`.
#'
#'This function is used for model comparison: does a 1- or 2-compartment TK
#'model fit the data any better than this naive "flat" model?
#'
#'(A note on nomenclature: Strictly speaking, in this model, concentration is
#'linear with dose. It is called "flat" because we usually consider
#'concentration normalized to dose, which is flat -- equal to the constant `A`.)
#'
#'@param params A named list of model parameter values. Must include: an item
#'  named `A`, representing the average ratio of concentration to dose.
#'@param time A numeric vector of times in hours.
#'@param dose A numeric vector of doses in mg/kg
#'@param iv.dose A logical vector: TRUE for single IV bolus dose; FALSE for
#'  single oral dose. Not used, but must be present for compatibility with other
#'  model functions.
#'
#'@return A vector of plasma concentration values (mg/L) corresponding to
#'  \code{time}.
#'
#'@author Caroline Ring, John Wambaugh, Chris Cook
#'
#'@export cp_flat

cp_flat <- function(time, params, dose, iv.dose) {

  cp <- params$A * dose

  return(cp)
}

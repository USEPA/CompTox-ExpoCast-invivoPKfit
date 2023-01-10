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
#'@param params A named list of model parameter values. Must include an item
#'  named `A`, representing the average ratio of concentration to dose. If any
#'  `medium %in% 'blood'`, then `params` must also include an item named
#'  `Rblood2plasma`, i.e. the ratio of chemical concentration in whole blood to
#'  the chemical concentration in plasma.
#'@param time A numeric vector of times in hours.
#'@param dose A numeric vector of doses in mg/kg
#'@param iv.dose A logical vector: TRUE for single IV bolus dose; FALSE for
#'  single oral dose. Not used, but must be present for compatibility with other
#'  model functions.
#'@param medium A character vector reflecting the medium in which each resulting
#'  concentration is to be calculated: "blood" or "plasma". Default is "plasma".
#'  Must be same length as other arguments, or length 1.
#'
#'@return A vector of plasma concentration values (mg/L) corresponding to
#'  \code{time}.
#'
#'@author Caroline Ring, John Wambaugh, Chris Cook
#'
#'@export cp_flat

cp_flat <- function(time, params, dose, iv.dose, medium) {

  cp <- params$A * dose
  cp[medium %in% "blood"] <- params$Rblood2plasma * cp[medium %in% "blood"]

  return(cp)
}

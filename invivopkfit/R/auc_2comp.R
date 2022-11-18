#' Analytical AUC for the 2-compartment model
#'
#' Calculate area under the plasma concentration vs. time curve for the
#' 1-compartment model, using an analytical equation (the integral of the
#' 1-compartment model equation with respect to time).
#'
#' @param params A named list of parameter values including the following:
#'   \describe{ \item{k12}{Rate at which compound moves from central to
#'   peripheral compartment} \item{k21}{Rate at which compound moves from
#'   peripheral to central compartment} \item{kelim}{Elimination rate}
#'   \item{V1}{Apparent volume of central compartment} } For oral administration
#'   (\code{iv.dose} FALSE), \code{params} must also include: \describe{
#'   \item{Fgutabs}{Oral bioavailability} \item{kgutabs}{Rate of absorption from
#'   gut} }
#'
#' @author Caroline Ring, John Wambaugh
#' @param time A vector of time values, in hours
#' @param dose A dose in mg/kg
#' @param iv.dose TRUE for single IV bolus dose, FALSE for single oral dose
#'@return A vector of plasma AUC values, evaluated at each time point in `time`.
#' @export auc_2comp
auc_2comp <- function(params, time, dose, iv.dose){

  if (any(sapply(params,function(x) identical(x,numeric(0))))) return(NA)

  if (is.null(params$Fgutabs) | is.na(params$Fgutabs))
  {
    params$Fgutabs <- 1
  }
  if (is.null(params$kgutabs) | is.na(params$kgutabs))
  {
    params$kgutabs <- 1
  }
  if (is.null(params$Fbetaofalpha)) params$Fbetaofalpha <- 1
  if (is.null(params$Ralphatokelim)) params$Ralphatokelim <- 1

  if (is.na(params$Fbetaofalpha)) params$Fbetaofalpha <- 1
  if (is.na(params$Ralphatokelim)) params$Ralphatokelim <- 1

  if (params$Fgutabs > 1) params$Fgutabs <- 1
  if (params$Fbetaofalpha > 1) params$Fbetaofalpha <- 1
  if (params$Ralphatokelim < 1) params$Ralphatokelim <- 1

  alpha <- params$Ralphatokelim * (params$kelim + 10^-6)
  beta <- params$Fbetaofalpha * alpha

  # try to keep k21 and k12 positive:
  k21 <- max(min(alpha * beta / params$kelim, alpha + beta - params$kelim), 0)
  k12 <- alpha + beta - params$kelim - k21

  alphabeta.sum <- alpha + beta
  alphabeta.prod <- alpha * beta

  if (iv.dose){ #for IV dosing
    A <- (dose * (alpha - k21)) / (params$V1 * (alpha - beta))
    B <- (dose * (k21 - beta))/(params$V1 * (alpha - beta))
    auc <- A/alpha * (1 - exp(-time*alpha)) +
      B/beta * (1 - exp(-time*beta))
  }else{
    A <- (params$Fgutabs * dose * (alpha - k21)) / (params$V1 * (alpha - beta))
    B <- (params$Fgutabs * dose * (k21 - beta)) / (params$V1 * (alpha - beta))

    auc <- A/alpha * (1 - exp(-time*alpha)) +
      B/beta * (1 - exp(-time*beta)) +
      (-A - B)/ka * (1 - exp(-time*ka))
  }

  auc[auc<0] <- 0
  auc[!is.finite(auc)] <- 0

  return(auc)
}

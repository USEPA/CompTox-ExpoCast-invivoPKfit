#' Analytical 2-compartment model
#'
#' @param params A named list of parameter values including the following:
#' \describe{ \item{k12}{Rate at which compound moves from central to
#' peripheral compartment} \item{k21}{Rate at which compound moves from
#' peripheral to central compartment} \item{kelim}{Elimination rate}
#' \item{V1}{Apparent volume of central compartment} } For oral administration
#' (\code{iv.dose} FALSE), \code{params} must also include: \describe{
#' \item{Fgutabs}{Oral bioavailability} \item{kgutabs}{Rate of absorption from
#' gut} }
#'
#' !author Caroline Ring, John Wambaugh
#' @param time A vector of time values, in hours
#' @param dose A dose in mg/kg
#' @param iv.dose TRUE for single IV bolus dose, FALSE for single oral dose
#' @return A vector of plasma concentration values corresponding to each value
#' in \code{time}
#' @export cp_2comp
cp_2comp <- function(params, time, dose, iv.dose)
{

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

  cp <- A * exp(-alpha * time) + B * exp(-beta * time)
  }else{ #for oral dosing
  A <- (params$Fgutabs * dose * (alpha - k21)) / (params$V1 * (alpha - beta))
  B <- (params$Fgutabs * dose * (k21 - beta)) / (params$V1 * (alpha - beta))
  C <- -(A + B)
  cp <-  A * exp(-alpha * time) + B * exp(-beta * time) + C * exp(-params$kgutabs * time)
  }

  # cp[cp < 10^-20] <- 10^-20
  cp[cp<0] <- 0
  cp[!is.finite(cp)] <- 0

  return(cp)
}

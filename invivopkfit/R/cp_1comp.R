#'Analytical 1-compartment model
#'
#' Calculates plasma concentrations according to the analytical solution for the  1-compartment model.
#'
#'@param time A vector of times in hours
#'@param params A named list of model parameter values. Must include:
#'\describe{
#'\item{kelim}{Elimination rate, 1/h}
#'\item{Vdist}{Volume of distribution, L/kg body weight}}
#'For oral administration (\code{iv.dose} FALSE), \code{params} must also include:
#'\describe{
#'\item{Fgutabs}{Oral bioavailability, unitless fraction}
#'\item{kgutabs}{Oral absorption rate, 1/h}}
#'@param dose Dose in mg/kg
#'@param iv.dose TRUE for single IV bolus dose; FALSE for single oral dose
#'
#'@return A vector of plasma concentration values corresponding to \code{time}.
#'
#'@author Caroline Ring, John Wambaugh
#'
#' @export cp_1comp
cp_1comp <- function(time, params, dose, iv.dose){

  if (any(sapply(params,function(x) identical(x,numeric(0))))) return(0)

  if (is.null(params$Fgutabs)|is.na(params$Fgutabs))
  {
    params$Fgutabs<-1
  }
  if (is.null(params$kgutabs)|is.na(params$kgutabs))
  {
    params$kgutabs<-1
  }
  if (params$Fgutabs>1) params$Fgutabs<-1

  if (iv.dose){

    cp <- dose*exp(-params$kelim * time)/params$Vdist

    }else{
      #params$Vdist <- log(params$Vdist) ##????? Why was this here and uncommented?
      const <- (params$Fgutabs * dose * params$kgutabs)/
        ((params$Vdist*(params$kgutabs - params$kelim)))

    cp <- const * (exp(-params$kelim * time) - exp(-params$kgutabs* time))

    #Note: this fails if kgutabs == kelim,
    #which usually happens if both are at the upper bound
  }
  cp[cp<0] <- 0
  cp[!is.finite(cp)] <- 0
  return(cp)
}

#' Analytical 1-compartment function
#'
#' Analytical 1-compartment function
#'
#'
#' @param params A named list of model parameter values. Must include:
#' \describe{ \item{kelim}{Elimination rate, 1/h} \item{Vdist}{Volume of
#' distribution, L/kg body weight}} For oral administration (\code{iv.dose}
#' FALSE), \code{params} must also include: \describe{ \item{Fgutabs}{Oral
#' bioavailability, unitless fraction} \item{kgutabs}{Oral absorption rate,
#' 1/h}}
#' @param dose Dose in mg/kg
#' @param times A vector of observation times in hours
#' @param time.units Either (h)ours or (d)ays (defaults to hours)
#' @param iv.dose TRUE for single IV bolus dose; FALSE for single oral dose
#' @return A matrix of \code{time}, plasma concentration, and area under the
#' curve (AUC).
#' @author Caroline Ring
#' @export analytic_1comp_fun
analytic_1comp_fun <- function(params, dose, times, time.units="h", iv.dose){

  #params: list with Vdist, Fgutabs, kgutabs, kelim
  #dose in units of mg/kg
  if (time.units=='d') times <- times*24 #convert times from days to hours
  #times in hours to match units of kgutabs and kelim (1/h)

  Cp <- cp_1comp(time=times,
                 params=params,
                 dose=dose,
                 iv.dose=iv.dose)

  #Cp will have units mg/L

  #get AUC
  #first try proper adaptive quadrature method
  AUC <-tryCatch(stats::integrate(cp_1comp,
                                  lower=min(times),
                                  upper=max(times),
                                  params=params,
                                  dose=dose,
                                  iv.dose=iv.dose)$value,
                 error = function(err){ #if integrate fails, use trapezoidal rule instead
                   cat('analytic_1comp_fun: Numerical integration failed; using trapezoidal rule\n')
                   return(sum(diff(times)*(Cp[-1]+Cp[-length(Cp)]))/2)
                 }
  )

  out.mat <- matrix(rep(0, length(times)*3),
                    nrow=length(times),
                    dimnames=list(NULL, c('time',
                                          'Ccompartment',
                                          'AUC')))
  if (time.units=='h') {
    out.mat[, 'time'] <- times
  } else if (time.units=='d'){
    out.mat[, 'time'] <- times/24 #convert from hours back to days
  }
  out.mat[, 'Ccompartment'] <- Cp
  out.mat[, 'AUC'] <- AUC

  return(out.mat)
}

#'Evaluate analytic models
#'
#'@param params A named list of parameter values. Must match the parameters of \code{model}.
#'@param dose A dose in units of mg/kg.
#'@param times A vector of times in hours or days. Ideally should be sorted.
#'@param time.units The units of \code{times}: "h" for hours, "d" for days.
#'@param iv.dose TRUE for IV dosing, FALSE for PO dosing.
#'@param model Analytic model to evaluate. Currently only "1compartment" or
#'  "2compartment" are implemented.
#'            
#' @author Caroline Ring
#'
#'@return A matrix with three columns:
#'\describe{
#'\item{time}{Time in days}
#'\item{Ccompartment}{Plasma concentrations}
#'\item{AUC}{Area under the curve}}
analytic_model_fun <- function(params,
                               dose,
                               times,
                               time.units,
                               iv.dose,
                               model){

  #params: list with Vdist, Fgutabs, kgutabs, kelim
  #dose in units of mg/kg
  if (time.units=='d') times <- times*24 #convert times from days to hours
  times <- sort(times) #Ensure times are sorted
  #times in hours to match units of kgutabs and kelim (1/h)

  #Assign model function to be evaluated
  mfun <- switch(model,
                 '1compartment' = cp_1comp,
                 '2compartment' = cp_2comp)

  Cp <- tryCatch(do.call(mfun,
                list(time=times,
              params=params,
              dose=dose,
              iv.dose=iv.dose)),
              error=rep(0,length=times))  

  #Cp will have units mg/L

  #Cumulative trapezoidal rule
  AUC <- c(0,cumsum(diff(times)*(Cp[-1]+Cp[-length(Cp)])))

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

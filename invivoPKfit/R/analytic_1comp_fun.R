#written by Caroline Ring
analytic_1comp_fun <- function(params, dose, times, time.units, iv.dose){
  
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
  AUC <-tryCatch(integrate(cp_1comp, 
                   lower=min(times), 
                   upper=max(times), 
                   params=params, 
                   dose=dose, 
                   iv.dose=iv.dose)$value,
    error = function(err){ #if integrate fails, use trapezoidal rule instead
      print('Numerical integration failed; using trapezoidal rule')
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
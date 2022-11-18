#' Compute model-predicted concentration values for a vector of time values
#'
#' Compute model-predicted concentration values for a vector of time values
#'
#'
#' @param design.times A vector of times from the design matrix
#' @param design.dose A dose from the design matrix
#' @param design.iv TRUE if IV data, FALSE if PO data
#' @param design.times.max The maximum of design.times
#' @param design.time.step The maximum number of time steps per hour.
#' @param modelfun "analytic" to use analytic model solution, "full" to use
#' full ODE model
#' @param model "1compartment" or "2compartment"
#' @param model.params A named vector of model parameter values.
#'
#' author Caroline Ring, John Wambaugh
#' @return A vector of the model-predicted plasma concentration values of the
#' same length as design.times, for each time in design.times
fitfun <- function(design.times,
                   design.dose,
                   design.iv,
                   design.times.max,
                   design.time.step,
                   modelfun,
                   model,
                   model.params){

  #Use the analytical solution
  if (modelfun == 'analytic'){
    #get sorted list of time points
    these.times <- sort(unique(design.times))
    try (out <- analytic_model_fun(params = model.params,
                                   times = these.times,
                                   time.units = 'd',
                                   dose=design.dose,
                                   iv.dose=design.iv,
                                   model=model))
    #Ways that this can fail:
    #It could fail with an error
    #Cp could be negative
    #Cp could be non-finite

    #if it fails with an error
    if (inherits(out, "try-error")){
      cat("fitfun: Error in analytic_model_fun\n")
      cat(paste(paste(apply(data.frame(Names=names(model.params),
                                       Values=model.params,
                                       stringsAsFactors=FALSE),
                            1,function(x) paste(x,collapse=": ")),collapse=", "),"\n",sep=""))
      browser()
      #      browser() #kick to debugger to find out what went wrong
      out[,"time"] <- these.times
      out[,"Ccompartment"] <- 10^-20 #just set concentration value to essentially zero
    }

    #if Cp is non-finite
    if (any(!is.finite(out[, 'Ccompartment']))){
      cat("fitfun: Error, Cp is non-finite\n")
      cat(paste(paste(apply(data.frame(Names=names(model.params),
                                       Values=model.params,
                                       stringsAsFactors=FALSE),
                            1,function(x) paste(x,collapse=": ")),collapse=", "),"\n",sep=""))
      #      browser() #kick to debugger to find out what went wrong
      out[,"time"] <- these.times
      out[,"Ccompartment"] <- 10^-20 #just set concentration value to essentially zero

    }

  } else if (modelfun == 'full'){
    #get sorted list of time points
    these.times <- sort(unique(c(0,
                                 design.times)))
    try(out <- httk::solve_1comp(parameters=model.params,
                           times=these.times,
                           dose=design.dose,
                           days=design.times.max,
                           tsteps = round(1/design.time.step),
                           output.units='mg/L',
                           initial.value=0,
                           iv.dose=design.iv,
                           suppress.messages=TRUE))

    if (inherits(out, "try-error")) #if solving 1-compartment model fails
    {
      out[,"time"] <- these.times
      out[,"Ccompartment"] <- 10^-20 #just set concentration value to essentially zero
      browser() #kick to debugger to find out what went wrong
    }

    #If concentration is negative, try again with lower tolerance
    while(any(out[,"Ccompartment"]<0) & atol > 1e-20)
    {
      atol <- atol/50
      try(out <- httk::solve_1comp(parameters=model.params,
                             times=these.times,
                             dose=design.dose,
                             days=design.times.max,
                             tsteps = round(1/design.time.step),
                             output.units='mg/L',
                             initial.value=0,
                             iv.dose=design.iv,
                             suppress.messages=TRUE,
                             atol=atol))
      if (inherits(out, "try-error"))
      {
        out[,"time"] <- these.times
        out[,"Ccompartment"] <- 10^-20
        browser()
      }
    }
  }


  #Result is "out": a matrix with three columns:
  #time, Ccompartment, and AUC
  #return a vector of the model-predicted Ccompartment values
  #of the same length as design.times, for each time in design.times

  #I think the easiest way to do that is to convert it into a data.table,
  #key by time column,
  #and then index by these.times
  out.dt <- as.data.table(out)
  setkey(out.dt, time)
  pred <- out.dt[J(design.times), Ccompartment]

  return(pred)
} #endfitfun

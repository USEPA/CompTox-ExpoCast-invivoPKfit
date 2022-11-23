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

    out <- tryCatch(analytic_model_fun(params = model.params,
                                   times = these.times,
                                   time.units = 'd',
                                   dose=design.dose,
                                   iv.dose=design.iv,
                                   model=model),
                     error = function(err){
                       cat("fitfun: Error in analytic_model_fun\n")
                       cat(paste0(err$message, "\n"))
                       cat(paste(paste(apply(data.frame(Names=names(model.params),
                                                        Values=model.params,
                                                        stringsAsFactors=FALSE),
                                             1,
                                             function(x) paste(x,collapse=": ")),
                                       collapse=", "),
                                 "\n",
                                 sep=""))
                       out_tmp <- cbind("time" = these.times,
                                        "Ccompartment" = rep(NA_real_, length(these.times)),
                                        "AUC" = rep(NA_real_, length(these.times)))
                       out_tmp
                     })

  } else if (modelfun == 'full'){
    #get sorted list of time points
    these.times <- sort(unique(c(0,
                                 design.times)))
    out <- tryCatch(httk::solve_1comp(parameters=model.params,
                           times=these.times,
                           dose=design.dose,
                           days=design.times.max,
                           tsteps = round(1/design.time.step),
                           output.units='mg/L',
                           initial.value=0,
                           iv.dose=design.iv,
                           suppress.messages=TRUE),
                    error = function(err){
                      cat("fitfun: Error in httk::solve_1comp()\n")
                      cat(paste0(err$message, "\n"))
                      cat(paste(paste(apply(data.frame(Names=names(model.params),
                                                       Values=model.params,
                                                       stringsAsFactors=FALSE),
                                            1,function(x) paste(x,collapse=": ")),collapse=", "),"\n",sep=""))
                      out_tmp <- cbind("time" = these.times,
                                       "Ccompartment" = rep(NA_real_, length(these.times)),
                                       "AUC" = rep(NA_real_, length(these.times)))
                      out_tmp
                    })

    #If concentration is negative, try again with lower tolerance
    while(any(out[,"Ccompartment"]<0) & atol > 1e-20)
    {
      atol <- atol/50
      out <- tryCatch(httk::solve_1comp(parameters=model.params,
                                        times=these.times,
                                        dose=design.dose,
                                        days=design.times.max,
                                        tsteps = round(1/design.time.step),
                                        output.units='mg/L',
                                        initial.value=0,
                                        iv.dose=design.iv,
                                        suppress.messages=TRUE,
                                        atol = atol),
                      error = function(err){
                        cat("fitfun: Error in httk::solve_1comp()\n")
                        cat(paste0(err$message, "\n"))
                        cat(paste(paste(apply(data.frame(Names=names(model.params),
                                                         Values=model.params,
                                                         stringsAsFactors=FALSE),
                                              1,function(x) paste(x,collapse=": ")),collapse=", "),"\n",sep=""))
                        out_tmp <- cbind("time" = these.times,
                                         "Ccompartment" = rep(NA_real_, length(these.times)),
                                         "AUC" = rep(NA_real_, length(these.times)))
                        out_tmp
                      })
    }
  }


  #Result is "out": a matrix with three columns:
  #time, Ccompartment, and AUC
  #return a vector of the model-predicted Ccompartment values
  #indexed by design.times

  times_order <- match(design.times, out[, "time"])
  pred <- out[times_order, "Ccompartment"]
  return(pred)
} #endfitfun

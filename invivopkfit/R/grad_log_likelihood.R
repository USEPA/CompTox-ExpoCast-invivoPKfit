#' Analytical gradient for log-likelihood function.
#'
#' @param opt_params A named vector of log-scaled parameter values. When this
#'   function is the objective function for a numerical optimizer, these are the
#'   parameters to be optimized.
#' @param const_params A named vector of additional log-scaled parameter values.
#'   When this function is the objective function for a numerical optimizer,
#'   these are additional model parameters whose value is to be held constant
#'   while the other parameters are optimized. Default NULL (meaning that all
#'   model parameters are supplied in `opt_params`). (If you are calling this
#'   function directly, you probably want to leave `const_param = NULL` and
#'   just supply all model parameters in `opt_params`.)
#' @param DF A `data.frame` of concentration-time data
#' @param modelfun "analytic" to use the analytic model solution, "full" to use
#'   the full ODE model
#' @param model The model to fit. Currently, only "flat", "1compartment" or
#'   "2compartment" models are implemented.
#' @param force_finite Logical: Whether to force the function to return a finite
#'   log-likelihood (e.g., as required by [optimx::optimx()] with method
#'   'L-BFGS-B'.) Default FALSE, allowing the function to return -Inf for
#'   infinitely-unlikely parameter combinations. When `force_finite == TRUE`,
#'   the function will replace -Inf with -999999.
#' @return The gradient of the log-likelihood function.
grad_log_likelihood <- function(opt_params,
                                const_params = NULL,
                                DF,
                                modelfun,
                                model,
                                LOQ_factor = 2,
                                force_finite = FALSE){
  #combine parameters to be optimized and held constant,
  #and convert into a list, since that is what model functions expect
  params <- as.list(c(opt_params, const_params))

  #Back-transform the log-transformed parameters onto the natural scale
  params <- lapply(params, exp)

  if(any(is.na(params))) return(-Inf)

  model <- model

  if (model != "flat") {
    #If oral fraction absorbed is >100% for some reason, return -inf
    if (params[["Fgutabs"]] > 1) return(-Inf)
  }

  #Extract parameters whose names do not match 'sigma'
  #(that is, all the actual model parameters)
  model.params <- params[!grepl(x = names(params),
                                pattern = "sigma")]

  atol <- 1e-13 #set a tolerance

  #get predicted plasma concentration vs. time for the current parameter
  #values, by dose and route

  DF <- as.data.table(DF)

  DF[, pred := fitfun(design.times = Time.Days,
                      design.dose = unique(Dose),
                      design.iv = unique(iv),
                      design.times.max = unique(Max.Time.Days),
                      design.time.step = unique(Time.Steps.PerHour),
                      modelfun = modelfun,
                      model = model,
                      model.params = model.params),
     by = .(Dose, Route)]

  #Match sigmas to references:
  #get vector of sigmas, named as "sigma_ref_ReferenceID" or just "sigma"
  sigmas <- unlist(params[grepl(x=names(params),
                                pattern = "sigma")])
  nref <- length(sigmas)
  if(nref > 1){
    #get the Reference ID for each sigma, based on its name
    refs_sigmas <- gsub(x = names(sigmas),
                        pattern = "sigma_ref_",
                        replacement = "")
    #match the Reference ID and assign each sigma to its corresponding reference
    DF[, sigma.ref:=sigmas[match(Reference,
                                 refs_sigmas,
                                 nomatch = 0)]
    ]
  }else{ #if only one reference, the parameter is just called "sigma"
    DF[, sigma.ref := sigmas]
  }


  #Compute log-normal log-likelihood gradient (LLG)

  #Pre-computation:
  #log-scale concentration
  DF[!is.na(Value), y := log(Value)] #detects
  DF[is.na(Value), y:=log(LOQ * LOQ_factor)] #nondetects
  #log-scale predicted mean
  #add 1e-12 in case predicted mean is 0
  DF[pred==0, pred := pred + 1e-12]
  DF[, mu:=log(pred)]
  #standardize log-scale Value by log-scale mean & sigma.ref
  DF[, z:=(y - mu)/sigma.ref]

  #For detects: gradient of PDF
  #wrt mu:
  DF[!is.na(Value),
     llg_mu := z/sigma.ref]
  #wrt sigma:
  DF[!is.na(Value),
     llg_sigma := (z^2 - 1)/sigma.ref]

  #For non-detects: gradient of CDF
  #wrt mu:
  DF[is.na(Value),
     llg_mu := -dnorm(z)/pnorm(z)]
  #wrt sigma:
  DF[is.na(Value),
     llg_sigma := -dnorm(z)*z/pnorm(z)]

  #Joint log-likelihood is sum of individual LLs....
  #Joint gradient is the same
  llg_mu <- DF[, sum(llg_mu)]
  llg_sigma <- DF[, sum(llg_sigma)]

  #return joint gradient as a vector
  return(c(llg_mu, llg_sigma))
}

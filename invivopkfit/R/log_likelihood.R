#' Log-likelihood
#'
#' The log-likelihood function (probability of data given model parameters).
#'
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
#'
#' @return A log-likelihood value for the data given the parameter values in
#'   params
log_likelihood <- function(opt_params,
                           const_params = NULL,
                           DF,
                           modelfun,
                           model,
                           LOQ_factor = 2,
                           force_finite = FALSE) {

  #combine parameters to be optimized and held constant,
  #and convert into a list, since that is what model functions expect
  params <- as.list(c(opt_params, const_params))

  #Back-transform the log-transformed parameters onto the natural scale
  params <- lapply(params, exp)

  # if(any(is.na(params))) return(-Inf)

  model <- model

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

  #if pred is exactly 0, add 1e-12
  #to avoid blowing up log transformation
  DF[pred==0, pred:= 1e-12]

  #Compute log-normal log-likelihood (LL):
  #For detects: PDF
  DF[!is.na(Value),
                 loglike := dnorm(x = log(Value),
                       mean = log(pred),
                       sd = sigma.ref,
                       log = TRUE)]
  #For non-detets: CDF
  DF[is.na(Value),
                 loglike := pnorm(q = log(LOQ * LOQ_factor),
                       mean = log(pred),
                       sd = sigma.ref,
                       log.p = TRUE)]

  #And sum over references to get overall LL
  ll <- DF[, sum(loglike)]
  #do *not* remove NAs, because they mean this parameter combination is impossible!

  #If ll isn't finite -- for example if a predicted concentration was negative --
  #just set it to -Inf to indicate that these parameters are infinitely unlikely
  if (!is.finite(ll)) ll <- -Inf

  #If user has selected to force return of a finite value,
  #e.g. as required by optimix with method 'L-BFGS-B',
  #then when log-likelihood is infinitely unlikely,
  #return a large negative number instead
  if(force_finite %in% TRUE){
  if (!is.finite(ll)) ll <- -999999
  }

  return(ll)
}

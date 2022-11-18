#' The log-likelihood function used for fitting
#'
#' The log-likelihood function used for fitting
#'
#'
#' @param params A named list of log-scaled parameter values
#' @param DF A `data.frame` of Concentration-time data for a given chemical
#' @param modelfun "analytic" to use the analytic model solution, "full" to use
#' the full ODE model
#' @param model The model to fit, either "1compartment" or "2compartment"
#' (other models not implemented for now)
#'
#' @return A log-likelihood value for the data given the parameter values in
#' params
log_likelihood <- function(params,
                           DF,
                           modelfun,
                           model) {

  #Back-transform the log-transformed parameters (see analyze_pk_data) into
  #their original form
  params <- lapply(params, exp)

  if(any(is.na(params))) return(-Inf)

  model <- model

  if (model != "flat") {
    #If oral fraction absorbed is >100% for some reason, set it to 100%
    if (params[["Fgutabs"]] > 1) return(-Inf)
  }

  #Extract parameters whose names do not match 'sigma'
  #(that is, all the actual model parameters)
  model.params <- params[!grepl(x = names(params),
                                pattern = "sigma")]

  atol <- 1e-13 #set a tolerance

  #get predicted plasma concentration vs. time for the current parameter
  #values, by dose and route
# if (nrow(DT) > 7) browser()
  # ### if the subset of data has LOQ values of 'NA', make the LOQ values = 0.45 * the minimum Value
  # DT[, LOQ := as.numeric(LOQ)]
  # DT[is.na(LOQ), LOQ := 0.45 * min(Value, na.rm = TRUE), by = .(Compound, Reference, Media)]
  # # , by = "Reference"]

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
  sigmas <- unlist(params[grepl(x=names(params),
                         pattern = "sigma")])
  nref <- length(sigmas)
  if(nref > 1){
    #get the Reference ID for each sigma, based on its name
    #sigma names are of pattern "sigma_ref_[Reference ID]", without the brackets
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

  #Compute log-normal log-likelihood (LL):
  # pred and sigma.ref are on the arithmetic scale, so convert to log-norma mu and var:

  DF[, mu := log((pred + 10 ^ -12) / (1 + sigma.ref ^ 2 / (pred + 10 ^ -12) ^ 2) ^ (1/2))]
  DF[, var := log(1 + sigma.ref ^ 2 / (pred + 10 ^ -12) ^ 2)]
  # Contibutions from detects
  ll.term1 <- DF[!is.na(Value),
                 sum(-((log(Value + 10 ^ -12) - mu) ^ 2 / 2 / var)
                     - log(Value*sqrt(2*var*pi)))]
  # Add in cumulative distribution up to twice LOQ for observations below 2*LOQ:
  ll.term2 <- DF[is.na(Value), sum(stats::plnorm( 2 * LOQ + 10 ^ -12, mean = mu, sd = var ^ (1/2), log.p = T))]

  #And sum over references to get overall LL
  ll <- ll.term1 + ll.term2

  #If ll isn't finite -- for example if a predicted concentration was negative --
  #just set it to -Inf to indicate that these parameters are infinitely unlikely
  if (!is.finite(ll)) ll <- -Inf


  # To improve estimation by effectively bounding the parameter estimates:
  # If any rates are for some reason absurdly high, reduce the likelihood:
  #  if (params[["kelim"]]>MAX.RATE) ll<-ll+1-params[["kelim"]]/MAX.RATE
  #  if (params[["kgutabs"]]>MAX.RATE) ll<-ll+1-params[["kgutabs"]]/MAX.RATE
  #  if ("k12" %in% names(params))
  #  {
  #    if (params[["k12"]]>MAX.RATE) ll<-ll+1-params[["k12"]]/MAX.RATE
  #    if (params[["k21"]]>MAX.RATE) ll<-ll+1-params[["k21"]]/MAX.RATE
  #  }
  return(ll)
}

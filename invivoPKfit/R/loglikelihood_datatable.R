#written by Caroline Ring
#modified by John Wambaugh
#
#' The log-likelihood function used for fitting
#'
#' @param params A named list of parameter values
#' @param DT Concentration-time data for a given chemical
#'@param modelfun "analytic" to use the analytic model solution, "full" to use
#'  the full ODE model
#'@param model The model to fit, either "1compartment" or "2compartment" (other
#'  models not implemented for now)
#'
#'  @return A log-likelihood value for the data given the parameter values in params

log.likelihood<- function(params,
                                DT,
                                modelfun,
                          model)
{
  # Bound the rates:
#  MAX.RATE <- 1000

  #Take a copy of the input data table so it behaves as though passed by value
  DT <- copy(DT)
  #Back-transform the log-transformed parameters (see analyze_pk_data) into
  #their original form
  params <- lapply(params,exp)

  if (any(unlist(lapply(params,is.na)))) return(-Inf)

  #If oral fraction absorbed is >100% for some reason, set it to 100%
  if (params[["Fgutabs"]]>1) return(-Inf)
  
#  #If a two-compartment model, then                alpha and beta must be real:
#  if ("k12" %in% names(params))
#  {
#    alphabeta.sum <- params[["kelim"]] + params[["k12"]] + params[["k21"]]
#    alphabeta.prod <- params[["kelim"]] * params[["k21"]]
#
#    alpha <- (alphabeta.sum + sqrt(alphabeta.sum^2 - 4*alphabeta.prod))/2
#    beta <- (alphabeta.sum - sqrt(alphabeta.sum^2 - 4*alphabeta.prod))/2
#    if (is.na(alpha) | is.na(beta)) return(-Inf)
#  }

  #Extract parameters whose names do not match 'sigma2'
  #(that is, all the actual model parameters)
  model.params <- params[names(params)[regexpr(pattern="sigma2",
                                               text=names(params))==-1]]

  DT[, sigma.ref:=params[paste("sigma2",
                                          Reference,
                                          sep=".")],
          by=Reference]

  atol <- 1e-13 #set a tolerance

    #get predicted plasma concentration vs. time for the current parameter
    #values, by dose and route
    DT[, pred:=fitfun(design.times = Time.Days,
                      design.dose = unique(Dose),
                      design.iv = unique(iv),
                      design.times.max = unique(Max.Time.Days),
                      design.time.step = unique(Time.Steps.PerHour),
                      modelfun = modelfun,
                      model=model,
                      model.params=model.params),
       by=.(Dose, Route)]

  #Construct the log-likelihood as Gaussian noise
  #assume residuals from each study are iid Gaussian

    #Compute LL terms by each Reference
    ll.term1 <- DT[!is.na(Value),
                       sum(-((log(Value+10^-12)-log(pred+10^-12))^2/
                             (2*sigma.ref)))]
    ll.term2 <-  DT[!is.na(Value)            , -(length(Value)*
          log(sqrt(sigma.ref*2*pi))),
    by=Reference]
    # Add in contributions from observations below 2*LOQ:
    ll.term3 <- DT[is.na(Value),sum(pnorm(log(2*LOQ+10^-12),mean=log(pred+10^-12),sd=sigma.ref^(1/2),log.p=T))]

    #And sum over references to get overall LL
    ll <- ll.term1 + ll.term2[, sum(V1)] + ll.term3

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

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


  #1-compartment model only....so far
  Fgutabs <- params$Fgutabs
  kgutabs <- params$kgutabs
  kelim <- params$kelim
  Vdist <- params$Vdist
  ##STILL need to add derivs wrt sigma_ref!!!

  #For detects
  #oral data
  #Detects, oral data, kelim
  DF[!is.na(Value) &
       Route %in% "po",
     grad_kelim := z*(Time*(kelim - kgutabs)*
                        exp(Time*(kelim + kgutabs)) +
                        (-exp(Time*kelim) +
                           exp(Time*kgutabs))*
                        exp(Time*kelim))*
       exp(-Time*kelim)/
       (sigma.ref*(kelim - kgutabs)*(exp(Time*kelim) - exp(Time*kgutabs)))
]
  #Detects, oral data, Vdist
  DF[!is.na(Value) &
       Route %in% "po",
     grad_Vdist := -z/(Vdist*sigma.ref)]
  #Detects, oral data, Fgutabs
  DF[!is.na(Value) &
       Route %in% "po",
     grad_Fgutabs := z/(Fgutabs*sigma.ref)]
  #Detects, oral data, kgutabs
  DF[!is.na(Value) &
       Route %in% "po",
     grad_kgutabs := z*(-Time*kelim*kgutabs*exp(Time*kelim) +
                          Time*kgutabs^2*exp(Time*kelim) +
                          kelim*exp(Time*kelim) -
                          kelim*exp(Time*kgutabs))/
       (kgutabs*sigma.ref*(kelim*exp(Time*kelim) -
                             kelim*exp(Time*kgutabs) -
                             kgutabs*exp(Time*kelim) +
                             kgutabs*exp(Time*kgutabs)))
     ]



  #Detects, iv data, kelim:
DF[!is.na(Value) & Route %in% "iv",
   grad_kelim := -Time*z/sigma.ref]

#Detects, IV data, Vdist:
DF[!is.na(Value) & Route %in% "iv",
   grad_Vdist:= -z/(Vdist*sigma.ref)]

#Detects, IV data, Fgutabs
DF[!is.na(Value) & Route %in% "iv",
   grad_Fgutabs:= 0]

#Detects, IV data, kgutabs
DF[!is.na(Value) & Route %in% "iv",
   grad_kgutabs:= 0]


#Nondetects, PO data, kelim
DF[is.na(Value) & Route %in% "po",
   grad_kelim := -dnorm(z)*(Time*(kelim - kgutabs)*
                         exp(Time*(kelim + kgutabs)) -
                         (exp(Time*kelim) - exp(Time*kgutabs))*
                         exp(Time*kelim))*exp(-Time*kelim)/
     (pnorm(z)*(kelim - kgutabs)*(exp(Time*kelim) - exp(Time*kgutabs)))]

#Nondetects, PO data, Vdist
DF[is.na(Value) & Route %in% "po",
   grad_Vdist := dnorm(z)/(pnorm(z)*Vdist)]

#Nondetects, PO data, Fgutabs:
DF[is.na(Value) & Route %in% "po",
   grad_Fgutabs := -dnorm(z)/(Fgutabs*pnorm(z))]

#Nondetcts, PO data, kgutabs:
DF[is.na(Value) & Route %in% "po",
   grad_kgutabs := dnorm(z)*(Time*kelim*kgutabs*exp(Time*kelim) -
                          Time*kgutabs^2*exp(Time*kelim) -
                          kelim*exp(Time*kelim) +
                          kelim*exp(Time*kgutabs))/
     (pnorm(z)*kgutabs*(kelim*exp(Time*kelim) -
                     kelim*exp(Time*kgutabs) -
                     kgutabs*exp(Time*kelim) +
                     kgutabs*exp(Time*kgutabs)))]

#Nondetects, IV data, kelim
DF[is.na(Value) & Route %in% "iv",
   grad_kelim := Time*dnorm(z)/pnorm(z)]

#Nondetects, IV data, Vdist:
DF[is.na(Value) & Route %in% "iv",
   grad_Vdist := dnorm(z)/(pnorm(z)*Vdist)]

#Nondetcts, IV data, Fgutabs
DF[is.na(Value) & Route %in% "iv",
   grad_Fgutabs := 0]

#Nondetcts, IV data, kgutabs
DF[is.na(Value) & Route %in% "iv",
   grad_kgutabs := 0]

### Sigma.ref ###

sigma_cols <- paste0("grad_",
                     grep(x = names(opt_params),
                          pattern = "sigma",
                          value = TRUE))

if(length(sigma_cols)==1){
  #if there is only one sigma,
  DF[!is.na(Value),
     (sigma_cols) := (z^2 - 1)/sigma.ref]

  DF[is.na(Value),
     (sigma_cols) := -dnorm(z)*z/pnorm(z)]
}else if(length(sigma_cols)>1){
  refs <- DF[, unique(Reference)]
  #assign for each reference
  for(this_ref in refs){
    this_sigma_col <- paste0("grad_sigma_ref_",
                             this_ref)
    other_sigma_cols <- setdiff(sigma_cols,
                                this_sigma_col)
    DF[Reference %in% this_ref &
         !is.na(Value),
       (this_sigma_col) := (z^2 - 1)/sigma.ref]

    DF[Reference %in% this_ref &
         is.na(Value),
       (this_sigma_col) := -dnorm(z)*z/pnorm(z)]

    #set the other sigma cols to zero
    DF[Reference %in% this_ref,
       (other_sigma_cols) := 0]
  }
}


ll <- DF[, sapply(.SD, sum),
         .SDcols= c("grad_kelim",
                    "grad_Vdist",
                    "grad_Fgutabs",
                    "grad_kgutabs",
                    grep(x = names(DF),
                         pattern = "grad_sigma",
                         value = TRUE))]

names(ll) <- gsub(x = names(ll),
                  pattern = "grad_",
                  replacement = "")


  #return gradient as a vector
  #in the same order as opt_params
  return(ll[names(opt_params)])
}

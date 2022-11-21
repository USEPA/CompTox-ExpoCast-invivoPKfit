hess_log_likelihood <- function(){
  function(opt_params,
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

    #for detects
    #second derivatives of log PDF
    #dmu, dmu:
    DF[!is.na(Value), hess11 := -1/sigma.ref^2]
    #dmu, dsigma:
    DF[!is.na(Value), hess12 :=  -2*z/sigma.ref^2]
    #dsigma, dmu
    DF[!is.na(Value), hess21 := -2*z/sigma.ref^2 ]
    #dsigma, dsigma
    DF[!is.na(Value), hess22 := (1 - 3*z^2)/sigma.ref^2 ]

    #for nondetects
    #second derivatives of log CDF
    #dmu, dmu:
    DF[is.na(Value), hess11 := -dnorm(z)*(pnorm(z)*z +
                                            dnorm(z)*sigma.ref)/
         (pnorm(z)^2*sigma.ref)
    ]
    #dmu, dsigma:
    DF[is.na(Value), hess12 :=  dnorm(z)*(-pnorm(z)*z^2 + pnorm(z) -
                                            dnorm(z)*sigma.ref*z)/
         (pnorm(z)^2*sigma.ref) ]
    #dsigma, dmu
    DF[is.na(Value), hess21 := dnorm(z)*(-pnorm(z)*z^2 + pnorm(z) -
                                           dnorm(z)*sigma.ref*z)/
         (pnorm(z)^2*sigma.ref)  ]
    #dsigma, dsigma
    DF[is.na(Value), hess22 := -dnorm(z)*z*(pnorm(z)*(z^2 - 2) +
                                              dnorm(z)*sigma.ref*z)/
         (pnorm(z)^2*sigma.ref) ]

    #sums
    hess_vals <- DF[, sapply(.SD, sum), .SDcols = c("hess11",
                                                    "hess12",
                                                    "hess21",
                                                    "hess22")]
    #form Hessian matrix
    hess_mat <- matrix(hess_vals,
                       nrow =2,
                       ncol =2,
                       byrow = TRUE)

    return(hess_mat)

    }

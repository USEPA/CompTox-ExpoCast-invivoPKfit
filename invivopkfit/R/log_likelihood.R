#' Log-likelihood
#'
#' The log-likelihood function (probability of data given model parameters).
#'
#'
#' @param params A named vector of log-scaled parameter values. When this
#'   function is the objective function for a numerical optimizer, these are the
#'   parameters to be optimized.
#' @param const_params A named vector of additional log-scaled parameter values.
#'   When this function is the objective function for a numerical optimizer,
#'   these are additional model parameters whose value is to be held constant
#'   while the other parameters are optimized. Default NULL (meaning that all
#'   model parameters are supplied in `params`). (If you are calling this
#'   function directly, you probably want to leave `const_params = NULL` and
#'   just supply all model parameters in `params`.)
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
log_likelihood <- function(params,
                           const_params = NULL,
                           DF,
                           modelfun,
                           model,
                           LOQ_factor = 2,
                           force_finite = FALSE) {

  #combine parameters to be optimized and held constant,
  #and convert into a list, since that is what model functions expect
  params <- as.list(c(params, const_params))

  model <- model

  #Extract parameters whose names do not match 'sigma'
  #(that is, all the actual model parameters)
  model.params <- params[!grepl(x = names(params),
                                pattern = "sigma")]

  atol <- 1e-13 #set a tolerance

  #get predicted plasma concentration vs. time for the current parameter
  #values, by dose and route

  if(modelfun %in% "analytic"){
  mfun <- switch(model,
                 '1compartment' = cp_1comp,
                 '2compartment' = cp_2comp,
                 'flat' = cp_flat)

  #get predicted plasma concentration
  pred <- do.call(mfun,
          args = list(params = model.params,
                      time = DF$Time, #in hours
                      dose = DF$Dose,
                      iv.dose = DF$iv
                      ))

  }else{ #evaluate full model
    if(!(model %in% "1compartment")){
      stop(paste0("modelfun = 'analytic' specfied with model = '",
                  model,
                  ". Full ODE model only implemented for 1-compartment model"))
    }else{
      #for each dose and route, eval httk::solve_1comp
      df_list <- split(DF, DF[c("Dose", "iv")])
      predlist <- lapply(df_list,
             function(this_df){
               these_times <- sort(unique(this_df$Time))
               this_tstep <- 1/min(diff(these_times))
               tmp <-  httk::solve_1comp(parameters=model.params,
                                          times=this_df$Time,
                                          dose=unique(this_df$Dose),
                                          days=max(this_df$Time.Days),
                                          tsteps = round(this_tstep),
                                          output.units = 'mg/L',
                                          initial.value=0,
                                          iv.dose = unique(this_df$iv),
                                          suppress.messages=TRUE)
               tmpdf <- data.frame(Ccomp = tmp[, 'Ccompartment'],
                                   Time = tmp[, "time"],
                                   Dose = unique(this_df$Dose),
                                   iv = unique(this_df$iv))
             })
      #bind list of data.frames into one big one
      pred_df <- Reduce(f = rbind, x = predlist)

      #add a tmp variable encoding row ordering in DF
      DF$rowid <- 1:.N

      #merge with DF
      pred_merge <- merge(DF, pred_df, by = c("Time", "Dose", "iv"))

      #sort rows in the same order as DF
      pred_merge <- pred_merge[order(pred_merge$rowid), ]

      #extract predicted plasma concentrations
      pred <- pred_merge$Ccomp

    }
  }

  #if pred is non-negative but tiny, then replace with machine tolerance
  pred[pred >= 0 &
         pred < .Machine$double.eps] <- .Machine$double.eps

  #Match sigmas to references:
  #get vector of sigmas, named as "sigma_ref_ReferenceID" or just "sigma"
  sigma_names <- grep(x = names(params),
                      pattern = "sigma",
                      fixed = TRUE,
                      value = TRUE)
  sigmas <- unlist(params[sigma_names])

  nref <- length(sigmas)
  if(nref > 1){
    #get the Reference ID for each sigma, based on its name
  refs_sigmas <- gsub(x = sigma_names,
                      pattern = "sigma_ref_",
                      replacement = "",
                      fixed = TRUE)
  #match the Reference ID and assign each sigma to its corresponding reference
  sigma_ref <- sigmas[match(DF$Reference,
                            refs_sigmas,
                            nomatch = 0)]

  }else{ #if only one reference, the parameter is just called "sigma"
    sigma_ref <- rep(sigmas, nrow(DF))
  }

  #log-transform predicted values -- suppress warnings about 'NaNs produced'
  suppressWarnings(logpred <- log(pred))

  #get log-likelihood for each observation
  loglike <- ifelse(is.na(DF$Value),
                    #for non-detects: CDF
                    pnorm(q = log(DF$LOQ * LOQ_factor),
                          mean = logpred,
                          sd = sigma_ref,
                          log.p = TRUE),
                    #for detects: PDF
                    dnorm(x = log(DF$Value),
                          mean = logpred,
                          sd = sigma_ref,
                          log = TRUE))
  #sum log-likelihoods over observations
  ll <- sum(loglike)
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

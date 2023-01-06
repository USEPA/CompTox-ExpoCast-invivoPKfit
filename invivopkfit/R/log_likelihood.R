#' Log-likelihood
#'
#' The log-likelihood function (probability of data given model parameters).
#'
#' The log-likelihood is formulated by assuming that residuals (model-predicted
#' concentrations minus observed concentrations) are independent and identically
#' distributed within a given combination of chemical, species, and reference,
#' obeying a zero-mean normal distribution with a chemical-, species-, and
#' reference-specific variance.
#'
#' # Log-likelihood equations
#'
#' For chemical-species combination \eqn{i} and reference \eqn{j}, define the
#' following quantities.
#'
#' \eqn{y_{ijk}} is the \eqn{k^{th}} observation of concentration, corresponding
#' to dose \eqn{d_{ijk}} and time \eqn{t_{ijk}}. Each observation has a
#' corresponding LOQ, \eqn{\textrm{LOQ}_{ijk}}.
#'
#' For multiple-subject observations, \eqn{y_{ijk}} is the \eqn{k^{th}} *sample
#' mean* observed concentration for chemical-species combination \eqn{i} and reference \eqn{j},
#' corresponding to dose \eqn{d_{ijk}} and time \eqn{t_{ijk}}. It represents the
#' mean of \eqn{n_{ijk}} individual measurements. It has a corresponding sample
#' standard deviation, \eqn{s_{ijk}}.
#'
#' \eqn{\bar{theta}_i} represents the vector of model parameters supplied in
#' argument `params` for chemical-species combination \eqn{i}.
#'
#' \eqn{\mu_{ijk}} is the model-predicted concentration for dose \eqn{d_{ijk}}
#' and time \eqn{t_{ijk}}. If \eqn{f(d, t; \bar{ \theta})} is the model function
#' evaluated at dose \eqn{d} and time \eqn{t}, with parameter vector
#' \eqn{\bar{\theta}}, then
#'
#'  \deqn{\mu_{ijk} = f \left(d_{ijk}, t_{ijk};
#' \bar{\theta}_i \right)}
#'
#' \eqn{\sigma_{ij}^2} is the reference- and chemical-specific residual
#' variance. (It is a hyperparameter.)
#'
#' ## Single-subject observations above limit of quantification (detects)
#'
#' This is the normal probability density function evaluated at the observed
#' concentration, as implemented in [stats::dnorm()].
#'
#' \deqn{LL_{ijk} = \log \left( \frac{1}{\sqrt{\sigma_{ij} 2 \pi}} \exp \left[
#' \frac{-1}{2} \left( \frac{y_{ijk} - \mu_{ijk}}{\sigma_{ij}} \right)^2
#' \right]   \right) }
#'
#'
#' ## Single-subject observations below limit of quantification (non-detects)
#'
#' This is the normal cumulative density function evaluated at the LOQ, as
#' implemented in [stats::pnorm()]. It is the total probability of observing a
#' concentration anywhere between 0 and the LOQ.
#'
#' \deqn{LL_{ijk} =  \log \left( \frac{1}{2} \left[ 1 + \textrm{erf} \left(
#' \frac{\textrm{LOQ}_{ijk} - \mu_{ijk}}{\sigma_{ij} \sqrt{2}} \right) \right]  \right) }
#'
#' ## Multiple-subject observations above limit of quantification
#'
#' This is the joint log-likelihood across the multiple subjects included in one
#' observation, re-expressed in terms of the sample mean, sample SD, and number
#' of subjects for that observation. It is implemented in [dnorm_summary()].
#'
#' \deqn{LL_{ijk} = n_{ijk} * \log \left[
#'  \frac{1}{\sigma_{ij} \sqrt{2 \pi}} +
#'   \frac{-1}{2 \sigma_{ij}^2}  \left(
#'    \left( n_{ijk} - 1 \right)
#'      s_{ijk}^2 + n_{ijk} \left( y_{ijk} - \mu_{ijk} \right)^2
#'       \right)
#'        \right]}
#'
#' ## Multiple-subject observations below limit of quantification
#'
#' This case is not implemented. If sample mean concentration is reported below
#' LOQ, then it is unclear what individual observed concentrations are
#' represented, and how they were combined to produce the summary data in terms
#' of sample mean and sample SD. Were all individual observations below LOQ? Or
#' were below-LOQ observations replaced with 0, LOQ/2, etc. before sample mean
#' and sample SD were computed? If the sample mean is reported below LOQ, what
#' LOQ is reported? Did individual observations all have the same LOQ, or is an
#' average or median LOQ being used?  It is impossible to formulate the
#' log-likelihood without knowing the answers to these questions. Therefore,
#' multiple-subject observations below LOQ are excluded from analysis (they are
#' removed in [preprocess_data()]).
#'
#' # Joint log-likelihood for a chemical and species
#'
#' The joint log-likelihood for a chemical and species is simply the sum of
#' log-likelihoods across observations and references.
#'
#' \deqn{LL_{i} = \sum_{j=1}^{J_i} \sum_{k=1}^{K_{ij}} LL_{ijk}}
#'
#' This is the overall probability of the observed data, given the model and
#' parameters.
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
#' @param DF A `data.frame` of concentration-time data, for example as produced
#'   by [preprocess_data()]. It must contain variables named `Reference`
#'   (representing the reference ID for each observation), `Time` (representing
#'   observation time points in hours), `Dose` (representing administered dose
#'   in mg/kg), `iv` (logical TRUE/FALSE, representing whether the administered
#'   dose was IV), `Value` (representing observed concentrations in mg/L), `LOQ`
#'   (numeric, representing limit of quantification of each observed
#'   concentration  in mg/L), and `N_Subjects` (representing the number of
#'   subjects included in the concentration measurement). If any `N_Subjects` >
#'   1, then `DF` must also include a variable `Value_SD`, representing the
#'   sample standard deviation for multiple-subject observations (where `Value`
#'   represents the sample mean of multiple measurements).
#' @param modelfun Character: whether to use the analytic TK model solution or
#'   solve the full ODE model numerically. "analytic" to use the analytic model
#'   solution, "full" to use the full ODE model. Default is "analytic".
#' @param model Character: The model to fit. Currently, only "flat",
#'   "1compartment" or "2compartment" models are implemented.
#' @param force_finite Logical: Whether to force the function to return a finite
#'   log-likelihood (e.g., as required by [optimx::optimx()] with method
#'   'L-BFGS-B'.) Default FALSE, allowing the function to return -Inf for
#'   infinitely-unlikely parameter combinations. When `force_finite == TRUE`,
#'   the function will replace -Inf with -999999.
#' @param negative Logical: Whether to return the *negative* log-likelihood
#'   (i.e., the log-likelihood multiplied by negative 1). Default TRUE, to
#'   return the negative log-likelihood. This option is useful when treating the
#'   log-likelihood as an objective function to be minimized by an optimization
#'   algorithm.
#'
#' @return A log-likelihood value for the data given the parameter values in
#'   params
log_likelihood <- function(params,
                           const_params = NULL,
                           DF,
                           modelfun = "analytic",
                           model,
                           force_finite = FALSE,
                           negative = TRUE) {

  #combine parameters to be optimized and held constant,
  #and convert into a list, since that is what model functions expect
  params <- as.list(c(params, const_params))

  #Extract parameters whose names do not match 'sigma'
  #(that is, all the actual model parameters)
  model.params <- params[!grepl(x = names(params),
                                pattern = "sigma")]

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
                                          days=ceiling(max(this_df$Time)/24),
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

  }else if(nref == 1){ #if only one reference, the parameter is just called "sigma"
    sigma_ref <- rep(sigmas, nrow(DF))
  }else{
    stop(paste("Could not find any parameters with 'sigma' in the name.",
    "Param names are:",
    paste(names(params), collapse = "; ")
    ))
  }


  #add sigma_ref and pred as temp columns to DF
  #this makes logical indexing easier
  DF$sigma_ref <- sigma_ref
  DF$pred <- pred

  #get log-likelihood for each observation

  #For single-subject observations:
  if(any(DF$N_Subjects %in% 1)){
  DF_single_subj <- subset(DF,
                           N_Subjects %in% 1)

  loglike_single_subj <- ifelse(is.na(DF_single_subj$Value),
                    #for non-detects: CDF
                    pnorm(q = DF_single_subj$LOQ,
                          mean = DF_single_subj$pred,
                          sd = DF_single_subj$sigma_ref,
                          log.p = TRUE),
                    #for detects: PDF
                    dnorm(x = DF_single_subj$Value,
                          mean = DF_single_subj$pred,
                          sd = DF_single_subj$sigma_ref,
                          log = TRUE))
  }else{
    loglike_single_subj <- 0
  }

  #For multi-subject observations:
  if(any(DF$N_Subjects > 1)){
  DF_multi_subj <- subset(DF, N_Subjects > 1)
  loglike_multi_subj <-  dnorm_summary(mu = DF_multi_subj$pred,
                                            sigma = DF_multi_subj$sigma_ref,
                                            x_mean = DF_multi_subj$Value,
                                            x_sd = DF_multi_subj$Value_SD,
                                            x_N = DF_multi_subj$N_Subjects,
                                            log = TRUE)
  }else{
    loglike_multi_subj <- 0
  }

  #sum log-likelihoods over observations
  ll <- sum(c(loglike_single_subj, loglike_multi_subj))
  #do *not* remove NAs, because they mean this parameter combination is impossible!

  #If ll isn't finite,
  #just set it to -Inf to indicate that these parameters are infinitely unlikely
  if (!is.finite(ll)) ll <- -Inf

  #if model predicted any negative concentrations,
  #then these parameters are infinitely unlikely
 if(any(pred < 0)) ll <- -Inf

  #If user has selected to force return of a finite value,
  #e.g. as required by optimx with method 'L-BFGS-B',
  #then when log-likelihood is infinitely unlikely,
  #return a large negative number instead
  if(force_finite %in% TRUE){
  if (!is.finite(ll)) ll <- -1 * 0.5*(.Machine$double.xmax)
  }

  #to get negative log-likelihood (e.g. for minimization)
  if(negative %in% TRUE){
    ll <- -1 * ll
  }

  return(ll)
}

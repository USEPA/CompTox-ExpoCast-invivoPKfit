#'Log-likelihood
#'
#'The log-likelihood function (probability of data given model parameters).
#'
#'The log-likelihood is formulated by assuming that residuals (model-predicted
#'concentrations minus observed concentrations) are independent and identically
#'distributed within a given combination of chemical, species, and study,
#'obeying a zero-mean normal distribution with a chemical-, species-, and
#'study-specific variance.
#'
#'# Log-likelihood equations
#'
#'For chemical-species combination \eqn{i} and study \eqn{j}, define the
#'following quantities.
#'
#'\eqn{y_{ijk}} is the \eqn{k^{th}} observation of concentration, corresponding
#'to dose \eqn{d_{ijk}} and time \eqn{t_{ijk}}. Each observation has a
#'corresponding LOQ, \eqn{\textrm{LOQ}_{ijk}}.
#'
#'For multiple-subject observations, \eqn{y_{ijk}} is the \eqn{k^{th}} *sample
#'mean* observed concentration for chemical-species combination \eqn{i} and
#'study \eqn{j}, corresponding to dose \eqn{d_{ijk}} and time \eqn{t_{ijk}}. It
#'represents the mean of \eqn{n_{ijk}} individual measurements. It has a
#'corresponding sample standard deviation, \eqn{s_{ijk}}.
#'
#'\eqn{\bar{theta}_i} represents the vector of model parameters supplied in
#'argument `params` for chemical-species combination \eqn{i}.
#'
#'\eqn{\mu_{ijk}} is the model-predicted concentration for dose \eqn{d_{ijk}}
#'and time \eqn{t_{ijk}}. If \eqn{f(d, t; \bar{ \theta})} is the model function
#'evaluated at dose \eqn{d} and time \eqn{t}, with parameter vector
#'\eqn{\bar{\theta}}, then
#'
#'  \deqn{\mu_{ijk} = f \left(d_{ijk}, t_{ijk};
#' \bar{\theta}_i \right)}
#'
#'\eqn{\sigma_{ij}^2} is the study- and chemical-specific residual variance. (It
#'is a hyperparameter.)
#'
#'## Single-subject observations above limit of quantification (detects)
#'
#'This is the normal probability density function evaluated at the observed
#'concentration, as implemented in [stats::dnorm()].
#'
#' \deqn{LL_{ijk} = \log \left( \frac{1}{\sqrt{\sigma_{ij} 2 \pi}} \exp \left[
#' \frac{-1}{2} \left( \frac{y_{ijk} - \mu_{ijk}}{\sigma_{ij}} \right)^2
#' \right]   \right) }
#'
#'
#'## Single-subject observations below limit of quantification (non-detects)
#'
#'This is the normal cumulative density function evaluated at the LOQ, as
#'implemented in [stats::pnorm()]. It is the total probability of observing a
#'concentration anywhere between 0 and the LOQ.
#'
#' \deqn{LL_{ijk} =  \log \left( \frac{1}{2} \left[ 1 + \textrm{erf} \left(
#' \frac{\textrm{LOQ}_{ijk} - \mu_{ijk}}{\sigma_{ij} \sqrt{2}} \right) \right]  \right) }
#'
#'## Multiple-subject observations above limit of quantification
#'
#'This is the joint log-likelihood across the multiple subjects included in one
#'observation, re-expressed in terms of the sample mean, sample SD, and number
#'of subjects for that observation. It is implemented in [dnorm_summary()].
#'
#' \deqn{LL_{ijk} = n_{ijk} * \log \left[
#'  \frac{1}{\sigma_{ij} \sqrt{2 \pi}} +
#'   \frac{-1}{2 \sigma_{ij}^2}  \left(
#'    \left( n_{ijk} - 1 \right)
#'      s_{ijk}^2 + n_{ijk} \left( y_{ijk} - \mu_{ijk} \right)^2
#'       \right)
#'        \right]}
#'
#'## Multiple-subject observations below limit of quantification
#'
#'This case is not implemented. If sample mean concentration is reported below
#'LOQ, then it is unclear what individual observed concentrations are
#'represented, and how they were combined to produce the summary data in terms
#'of sample mean and sample SD. Were all individual observations below LOQ? Or
#'were below-LOQ observations replaced with 0, LOQ/2, etc. before sample mean
#'and sample SD were computed? If the sample mean is reported below LOQ, what
#'LOQ is reported? Did individual observations all have the same LOQ, or is an
#'average or median LOQ being used?  It is impossible to formulate the
#'log-likelihood without knowing the answers to these questions. Therefore,
#'multiple-subject observations below LOQ are excluded from analysis (they are
#'removed in [preprocess_data()]).
#'
#'# Joint log-likelihood for a chemical and species
#'
#'The joint log-likelihood for a chemical and species is simply the sum of
#'log-likelihoods across observations and studys.
#'
#'\deqn{LL_{i} = \sum_{j=1}^{J_i} \sum_{k=1}^{K_{ij}} LL_{ijk}}
#'
#'This is the overall probability of the observed data, given the model and
#'parameters.
#'
#'@param params A named vector of log-scaled parameter values. When this
#'  function is the objective function for a numerical optimizer, these are the
#'  parameters to be optimized.
#'@param const_params A named vector of additional log-scaled parameter values.
#'  When this function is the objective function for a numerical optimizer,
#'  these are additional model parameters whose value is to be held constant
#'  while the other parameters are optimized. Default NULL (meaning that all
#'  model parameters are supplied in `params`). (If you are calling this
#'  function directly, you probably want to leave `const_params = NULL` and just
#'  supply all model parameters in `params`.)
#'@param DF A `data.frame` of concentration-time data, for example as produced
#'  by [preprocess_data()]. It must contain variables named `Study`
#'  (representing the study ID for each observation), `Time` (representing
#'  observation time points in hours), `Dose` (representing administered dose in
#'  mg/kg), `iv` (logical TRUE/FALSE, representing whether the administered dose
#'  was IV), `Value` (representing observed concentrations in mg/L), `LOQ`
#'  (numeric, representing limit of quantification of each observed
#'  concentration  in mg/L), and `N_Subjects` (representing the number of
#'  subjects included in the concentration measurement). If any `N_Subjects` >
#'  1, then `DF` must also include a variable `Value_SD`, representing the
#'  sample standard deviation for multiple-subject observations (where `Value`
#'  represents the sample mean of multiple measurements).
#'@param modelfun Character: whether to use the analytic TK model solution or
#'  solve the full ODE model numerically. "analytic" to use the analytic model
#'  solution, "full" to use the full ODE model. Default is "analytic".
#'@param model Character: The model to fit. Currently, only "flat",
#'  "1compartment" or "2compartment" models are implemented.]
#'@param fit_conc_dose Logical: Whether to fit dose-normalized concentrations
#'  (TRUE) or non-dose-normalized concentrations (FALSE). Default TRUE.
#'@param fit_log_conc TRUE or FALSE (default FALSE): whether to apply log-scaling to
#'  the concentrations (or dose-normalized concentrations) before evaluating the
#'  log-likelihood.
#'@param force_finite Logical: Whether to force the function to return a finite
#'  log-likelihood (e.g., as required by [optimx::optimx()] with method
#'  'L-BFGS-B'.) Default FALSE, allowing the function to return -Inf for
#'  infinitely-unlikely parameter combinations. When `force_finite == TRUE`, the
#'  function will replace -Inf with -999999.
#'@param negative Logical: Whether to return the *negative* log-likelihood
#'  (i.e., the log-likelihood multiplied by negative 1). Default TRUE, to return
#'  the negative log-likelihood. This option is useful when treating the
#'  log-likelihood as an objective function to be minimized by an optimization
#'  algorithm.
#'
#'@return A log-likelihood value for the data given the parameter values in
#'  params
log_likelihood <- function(params,
                           const_params = NULL,
                           DF,
                          mfun,
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

  pred <- do.call(mfun,
          args = list(params = model.params,
                      time = DF$Time, #in hours
                      dose = DF$Dose,
                      iv.dose = DF$iv,
                      medium = DF$Media
                      ))

  #Match sigmas to studys:
  #get vector of sigmas, named as "sigma_study_StudyID" or just "sigma"
  sigma_names <- grep(x = names(params),
                      pattern = "sigma",
                      fixed = TRUE,
                      value = TRUE)
  sigmas <- unlist(params[sigma_names])

  nstudy <- length(sigmas)
  if(nstudy > 1){
    #get the Study ID for each sigma, based on its name
  studies_sigmas <- gsub(x = sigma_names,
                      pattern = "sigma_study_",
                      replacement = "",
                      fixed = TRUE)
  #match the study ID and assign each sigma to its corresponding study
  sigma_study <- sigmas[match(DF$Study,
                            studies_sigmas,
                            nomatch = 0)]

  }else if(nstudy == 1){ #if only one study, the parameter is just called "sigma"
    sigma_study <- rep(sigmas, nrow(DF))
  }else{
    stop(paste("Could not find any parameters with 'sigma' in the name.",
    "Param names are:",
    paste(names(params), collapse = "; ")
    ))
  }


  #add sigma_study and pred as temp columns to DF
  #this makes logical indexing easier
  DF$sigma_study <- sigma_study
  DF$pred <- pred
  DF$pred_dose <- DF$pred/DF$Dose

  #if fit_conc_dose is TRUE, remove any zero-dose cases,
  #as these will have NaN log-likelihoods
  if(fit_conc_dose %in% TRUE){
    DF <- subset(DF, Dose > 0)
  }

  if(fit_conc_dose %in% TRUE){
    value <- DF$Value_Dose
    loq <- DF$LOQ_Dose
    pred <- DF$pred_dose
    value_sd <- DF$Value_SD_Dose
  }else{
    value <- DF$Value
    loq <- DF$LOQ
    pred <- DF$pred
    value_sd <- DF$Value_SD
  }

  if(fit_log_conc %in% TRUE){
    value_natural <- value
    value <- log(value)
    loq <- log(loq)
    pred_natural <- pred
    pred <- log(pred)
    ll_summary <- "dlnorm_summary"
  }else{
    ll_summary <- "dnorm_summary"
    value_natural <- value
    pred_natural <- pred
  }

  #get log-likelihood for each observation
  loglike <- ifelse(DF$N_Subjects %in% 1,
                    ifelse(DF$Detect %in% TRUE,
                           dnorm(x = value,
                                 mean = pred,
                                 sd = DF$sigma_study,
                                 log = TRUE),
                           pnorm(q = loq,
                                 mean = pred,
                                 sd = DF$sigma_study,
                                 log.p = TRUE)
                    ),
                    do.call(ll_summary,
                            list(mu = pred,
                                 sigma = DF$sigma_study,
                                 x_mean = value_natural,
                                 x_sd = value_sd,
                                 x_N = DF$N_Subjects,
                                 log = TRUE))
  )

  #sum log-likelihoods over observations
  ll <- sum(loglike)
  #do *not* remove NAs, because they mean this parameter combination is impossible!

  #If ll isn't finite,
  #just set it to -Inf to indicate that these parameters are infinitely unlikely
  if (!is.finite(ll)) ll <- -Inf

  #if model predicted any NA or Inf concentrations,
  #then parameters are infinitely unlikely
  if(any(!is.finite(pred_natural))) ll <- -Inf

  #if model predicted any negative concentrations,
  #then these parameters are infinitely unlikely
 if(any((pred_natural < 0) %in% TRUE)) ll <- -Inf

  #If user has selected to force return of a finite value,
  #e.g. as required by optimx with method 'L-BFGS-B',
  #then when log-likelihood is infinitely unlikely,
  #return a large negative number instead
  if(force_finite %in% TRUE){
  if (!is.finite(ll)) ll <- -1 * (.Machine$double.xmax)
  }

  #to get negative log-likelihood (e.g. for minimization)
  if(negative %in% TRUE){
    ll <- -1 * ll
  }

  return(ll)
}

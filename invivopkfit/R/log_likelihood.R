#'Log-likelihood
#'
#'The log-likelihood function (probability of data given model parameters).
#'
#'The log-likelihood is formulated by assuming that residuals (transformed
#'model-predicted concentrations minus transformed observed concentrations) are
#'independent, and that groups of residuals obey zero-mean normal distributions
#'where each group may have a separate error variance. Error groups are defined
#'as unique combinations of variables in the harmonized data, by a command such
#'as `pk(data = ...) + stat_error_model(error_group = vars(...)`.
#'
#'# Log-likelihood equations
#'
#'For chemical-species combination \eqn{i} and study \eqn{j}, define the
#'following quantities.
#'
#'\eqn{y_{ijk}} is the \eqn{k^{th}} observation of concentration (which may be
#'transformed, e.g. dose-normalized and/or log10-transformed), corresponding to
#'dose \eqn{d_{ijk}} and time \eqn{t_{ijk}}. Each observation has a
#'corresponding LOQ, \eqn{\textrm{LOQ}_{ijk}}.
#'
#'For multiple-subject observations, \eqn{y_{ijk}} is the \eqn{k^{th}} *sample
#'mean* observed concentration for chemical-species combination \eqn{i} and
#'study \eqn{j}, corresponding to dose \eqn{d_{ijk}} and time \eqn{t_{ijk}}. It
#'represents the mean of \eqn{n_{ijk}} individual measurements. It has a
#'corresponding sample standard deviation, \eqn{s_{ijk}}. In the harmonized
#'data,  \eqn{s_{ijk}} is contained in variable `Conc_SD`.
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
#'@param par A named list of parameters and their values that are being
#'  optimized.
#'@param const_params Optional: A named list of parameters and their values that
#'  are being held constant.
#'@param fitdata A `data.frame` of data with harmonized variable names, with transformations applied.
#'  Required: `Time_trans`, `Dose`, `Conc`, `Detect`, `N_Subjects`, `Conc_SD`,
#'  `Conc_trans`, `Conc_SD_trans`.
#'@param data_sigma_group A `factor` vector which could be a new variable in
#'  `fitdata`, giving the error group for each row in `fitdata`.
#'@param modelfun Character or function: The name of the function that produces
#'  concentration predictions for the model being evaluated.
#'@param scales_conc As from a [pk] object, element `$scales$conc` (see
#'  [scale_conc()] for how to specify. A list specifying the
#'  scaling/transformation of concentrations, with elements named
#'  `ratio_conc_dose`, `dose_norm`, `log10_trans`, `expr`.
#'@param force_finite Logical: Whether to force return of a finite value (e.g.
#'  as required by method `L-BFGS-B` in [optimx::optimx()]). Default FALSE. If
#'  TRUE, then if the log-likelihood works out to be non-finite, then it will be
#'  replaced with `.Machine$double.xmax`.
#'@param negative Logical: Whether to return the *negative* log-likelihood
#'  (i.e., the log-likelihood multiplied by negative 1). Default TRUE, to
#'  multiply the log-likelihood by negative 1 before returning it. This option
#'  is useful when treating the log-likelihood as an objective function to be
#'  *minimized* by an optimization algorithm.
#'@return A log-likelihood value for the data given the parameter values in
#'  params
#' @export
#' @author Caroline Ring
log_likelihood <- function(par,
                           const_params = NULL,
                           fitdata = NULL,
                           data_sigma_group = NULL,
                           modelfun = NULL,
                           scales_conc = NULL,
                           negative = TRUE,
                           force_finite = FALSE) {

  #combine parameters to be optimized and held constant,
  #and convert into a list, since that is what model functions expect
  params <- as.list(c(par, const_params))

  #Extract parameters whose names do not match 'sigma'
  #(that is, all the actual model parameters)
  model.params <- params[!grepl(x = names(params),
                                pattern = "sigma")]


  #get un-transformed predicted plasma concentration vs. time for the current parameter
  #values, by dose and route

  pred <- do.call(modelfun,
          args = list(params = model.params,
                      time = fitdata$Time_trans,
                      dose = fitdata$Dose,
                      route = fitdata$Route,
                      medium = fitdata$Media
                      ))

if(any(grepl(x = names(params),
             pattern= "sigma"))){
  sigma.params <- params[grepl(x = names(params),
                               pattern = "sigma")]
  #match the study ID and assign each sigma to its corresponding study
  sigma_study <- unlist(sigma.params[paste0("sigma_", data_sigma_group)])
  }else{
    stop(paste("Could not find any parameters with 'sigma' in the name.",
    "Param names are:",
    paste(names(params), collapse = "; ")
    ))
  }


  #add sigma_study and pred as temp columns to DF
  #this makes logical indexing easier
  fitdata$sigma_study <- sigma_study
  fitdata$pred <- pred
  #transform predictions the same way as concentrations
  fitdata$pred_trans <- rlang::eval_tidy(scales_conc$expr,
                                         data = cbind(fitdata,
                                                      data.frame(".conc" = pred))
  )


  if(scales_conc$log10_trans %in% TRUE){
    ll_summary <- "dlnorm_summary"
    #maintain other transformations such as dose-scaling,
    #but undo log10 transformation
    conc_natural <- 10^(fitdata$Conc_trans)
    conc_sd_natural <- 10^(fitdata$Conc_SD_trans)
  }else{ #if log10 transformation has *not* been applied
    ll_summary <- "dnorm_summary"
    conc_natural <- fitdata$Conc_trans
    conc_sd_natural <- fitdata$Conc_SD_trans
  }

  #get log-likelihood for each observation
  loglike <- ifelse(fitdata$N_Subjects %in% 1,
                    ifelse(fitdata$Detect %in% TRUE,
                           dnorm(x = fitdata$Conc_trans,
                                 mean = fitdata$pred_trans,
                                 sd = fitdata$sigma_study,
                                 log = TRUE),
                           pnorm(q = fitdata$Conc_trans,
                                 mean = fitdata$pred_trans,
                                 sd = fitdata$sigma_study,
                                 log.p = TRUE)
                    ),
                    do.call(ll_summary,
                            list(mu = fitdata$pred_trans,
                                 sigma = fitdata$sigma_study,
                                 x_mean = conc_natural,
                                 x_sd = conc_sd_natural,
                                 x_N = fitdata$N_Subjects,
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
  if(any(!is.finite(pred))) ll <- -Inf

  #if model predicted any negative concentrations,
  #then these parameters are infinitely unlikely
 if(any((pred < 0) %in% TRUE)) ll <- -Inf

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

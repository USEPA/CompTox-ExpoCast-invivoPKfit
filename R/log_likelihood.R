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
#'marked as excluded in [preprocess_data()]).
#'
#'# Joint log-likelihood for a chemical and species
#'
#'The joint log-likelihood for a chemical and species is simply the sum of
#'log-likelihoods across observations.
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
#'@param data A `data.frame` of data with harmonized variable names. Required:
#'  `Time_trans`, `Dose`, `Conc`, `Detect`, `N_Subjects`, `Conc_SD`. `Conc` and
#'  `Conc_SD` will be transformed according to `dose_norm` and `log10_trans`.
#'@param data_sigma_group A `factor` vector which could be a new variable in
#'  `data`, giving the error group for each row in `data`.
#'@param modelfun Character or function: The name of the function that produces
#'  concentration predictions for the model being evaluated.
#'@param dose_norm Logical: Whether to dose-normalize predicted and observed
#'  concentrations before calculating likelihood.
#'@param log10_trans Logical: Whether to apply a [log10()] transformation to
#'  predicted and observed concentrations before calculating likelihood.
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
#'@export
#'@author Caroline Ring
log_likelihood <- function(par,
                           const_params = NULL,
                           data = NULL,
                           data_sigma_group = NULL,
                           modelfun = NULL,
                           dose_norm = FALSE,
                           log10_trans = FALSE,
                           negative = TRUE,
                           force_finite = FALSE,
                           suppress.messages = TRUE) {
  #combine parameters to be optimized and held constant
  params <- c(par, const_params)

  #Extract parameters whose names do not match 'sigma'
  #(that is, all the actual model parameters)
  model.params <- params[!grepl(x = names(params),
                                pattern = "sigma")]

  #get un-transformed predicted plasma concentration vs. time for the current parameter
  #values, by dose and route
  pred <- do.call(
    modelfun,
    args = list(
      params = model.params,
      time = data$Time_trans,
      dose = data$Dose,
      route = data$Route,
      medium = data$Media
    )
  )


  data$pred <- pred

  #transform predictions and concentrations

  #dose normalization?
  if (dose_norm %in% TRUE) {
    data$pred_trans <- pred / data$Dose
    data$Conc_trans <- data$Conc / data$Dose
    data$Conc_SD_trans <- data$Conc_SD / data$Dose
  } else{
    data$pred_trans <- pred
    data$Conc_trans <- data$Conc
    data$Conc_SD_trans <- data$Conc_SD
  }

  #log10 transformation?
  if (log10_trans %in% TRUE) {
    data$pred_trans <- log10(data$pred_trans)
    data$Conc_trans <- log10(data$Conc_trans)
    data$Conc_SD_trans <- log10(data$Conc_SD_trans)
  } else{
    data$pred_trans <-  data$pred_trans
    data$Conc_trans <- data$Conc_trans
    data$Conc_SD_trans <- data$Conc_SD_trans
  }

  if (log10_trans %in% TRUE) {
    ll_summary <- "dlnorm_summary"
    #maintain other transformations such as dose-scaling,
    #but undo log10 transformation
    data$conc_natural <- 10 ^ (data$Conc_trans)
    data$conc_sd_natural <- 10 ^ (data$Conc_SD_trans)
  } else{
    #if log10 transformation has *not* been applied
    ll_summary <- "dnorm_summary"
    data$conc_natural <- data$Conc_trans
    data$conc_sd_natural <- data$Conc_SD_trans
  }

  #residual error SDs
  #defined by data_sigma_group
  if (any(grepl(x = names(params),
                pattern = "sigma"))) {
    #get sigma params
    sigma_params <- as.list(params[grepl(x = names(params),
                                         pattern = "sigma")])
    #match the study ID and assign each sigma to its corresponding study
    #except if data_sigma_group is NA -- then assume data are equally likely to come from any of the existing distributions
    sigma_obs <- sigma_params[paste("sigma",
                                    data_sigma_group,
                                    sep = "_")]
    sigma_obs <- sapply(sigma_obs,
                        function(x)
                          ifelse(is.null(x), NA_real_, x))
    data$sigma_obs <- sigma_obs



    if (any(!is.na(sigma_obs))) {
      #compute log likelihoods for observations with sigmas
      data_sigma <- data[!is.na(sigma_obs),]
      ll_data_sigma <-  ifelse(
        data_sigma$N_Subjects %in% 1,
        ifelse(
          data_sigma$Detect %in% TRUE,
          dnorm(
            x = data_sigma$Conc_trans,
            mean = data_sigma$pred_trans,
            sd = data_sigma$sigma_obs,
            log = TRUE
          ),
          pnorm(
            q = data_sigma$Conc_trans,
            mean = data_sigma$pred_trans,
            sd = data_sigma$sigma_obs,
            log.p = TRUE
          )
        ),
        do.call(
          ll_summary,
          list(
            mu = data_sigma$pred_trans,
            sigma = data_sigma$sigma_obs,
            x_mean = data_sigma$conc_natural,
            x_sd = data_sigma$conc_sd_natural,
            x_N = data_sigma$N_Subjects,
            log = TRUE
          )
        )
      )
    } else{
      ll_data_sigma <- numeric(0)
    }

    #compute log likelihoods for observations without sigmas
    if (any(is.na(sigma_obs))) {
      data_no_sigma <- data[is.na(sigma_obs),]
      if (suppress.messages %in% FALSE) {
        message(
          paste(
            "log_likelihood():",
            nrow(data_no_sigma),
            "observations are not in any existing error-SD (sigma) group.",
            "They will be treated as equally likely to be in",
            "any of the existing error-SD groups."
          )
        )
      }
      ll_data_no_sigma <-  sapply(unlist(sigma_params),
                                  function(this_sigma) {
                                    #place these observations in each group in turn
                                    ifelse(
                                      data_no_sigma$N_Subjects %in% 1,
                                      ifelse(
                                        data_no_sigma$Detect %in% TRUE,
                                        dnorm(
                                          x = data_no_sigma$Conc_trans,
                                          mean = data_no_sigma$pred_trans,
                                          sd = this_sigma,
                                          log = TRUE
                                        ),
                                        pnorm(
                                          q = data_no_sigma$Conc_trans,
                                          mean = data_no_sigma$pred_trans,
                                          sd = this_sigma,
                                          log.p = TRUE
                                        )
                                      ),
                                      do.call(
                                        ll_summary,
                                        list(
                                          mu = data_no_sigma$pred_trans,
                                          sigma = this_sigma,
                                          x_mean = data_no_sigma$conc_natural,
                                          x_sd = data_no_sigma$conc_sd_natural,
                                          x_N = data_no_sigma$N_Subjects,
                                          log = TRUE
                                        )
                                      )
                                    )
                                  }) #end sapply(unlist(sigma_params),

      #the result will be a matrix with as many columns as there are sigma_params,
      #and as many rows as there are observations without sigma groups
      #take the row means
      ll_data_no_sigma <- rowMeans(ll_data_no_sigma)
    } else{
      ll_data_no_sigma <- numeric(0)
    }

    loglike <- c(ll_data_sigma,
                 ll_data_no_sigma)

  } else{
    stop(
      paste(
        "Could not find any parameters with 'sigma' in the name.",
        "Param names are:",
        paste(names(params), collapse = "; ")
      )
    )
  }

  #sum log-likelihoods over observations
  ll <- sum(loglike)
  #do *not* remove NAs, because they mean this parameter combination is impossible!

  #If ll isn't finite,
  #just set it to -Inf to indicate that these parameters are infinitely unlikely
  if (!is.finite(ll))
    ll <- -Inf

  #if model predicted any NA or Inf concentrations,
  #then parameters are infinitely unlikely
  if (any(!is.finite(pred)))
    ll <- -Inf

  #if model predicted any negative concentrations,
  #then these parameters are infinitely unlikely
  if (any((pred < 0) %in% TRUE))
    ll <- -Inf

  #If user has selected to force return of a finite value,
  #e.g. as required by optimx with method 'L-BFGS-B',
  #then when log-likelihood is infinitely unlikely,
  #return a large negative number instead
  if (force_finite %in% TRUE) {
    if (!is.finite(ll)) {
      #now return sqrt of .Machine$double.xmax, not just .Machine$double.xmax
      #not taking sqrt seems to break L-BFGS-B sometimes
      ll <- -1 * sqrt(.Machine$double.xmax)
    }
  }

  #to get negative log-likelihood (e.g. for minimization)
  if (negative %in% TRUE) {
    ll <- -1 * ll
  }

  return(ll)
}

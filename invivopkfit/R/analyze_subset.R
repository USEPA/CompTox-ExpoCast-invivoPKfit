#' Fit a PK model to one set of concentration vs. time data
#'
#' Fit a specified PK model to one set of concentration vs. time data and return
#' a set of fitted parameter values.
#'
#' Typically this function is not called directly by the user, but is called
#' from [fit_all], the main fitting function.
#'
#' This function estimates one set of parameters for a specified PK model, based
#' on a set of concentration vs. time data in `fitdata`. `fitdata` is typically
#' data for one chemical and one species, but it may include data from more than
#' one study and more than one dosing route.
#'
#' The parameters to be estimated are those defined in
#' [get_opt_params()] for the model specified in `model` and `modelfun`:
#'
#'  - for `model = 'flat'`, the parameter estimated is `A`
#'  - for `model = '1compartment'`, the parameters estimated are `Vdist`,
#'  `kelim`, and possibly `kgutabs` and `Fgutabs` (see below).
#'  - for `model = '2compartment'`, the parameters estimated are `V1`,
#'  `Ralphatobeta`, `Fbetaofalpha`, and possibly `kgutabs` and `Fgutabs`.
#'
#'  For 1-compartment and 2-compartment models, `kgutabs` will be estimated from
#'  the data only if `fitdata` includes oral dosing data; otherwise it will be
#'  set to NA.
#'
#'  `Fgutabs` will be estimated from the data only if `fitdata` includes both
#'  oral and IV data. If `fitdata` includes oral data but not IV data, `Fgutabs`
#'  will be held constant at 1 while other parameters are estimated. If
#'  `fitdata` does not include oral data, `Fgutabs` will be set to NA.
#'
#' In addition to the model parameters, the log-scale standard deviation of the
#' residual errors will be estimated. If `fitdata` includes more than one unique
#' reference value (`length(unique(fitdata$Reference))>1`) and `pool_sigma ==
#' FALSE`, then a separate log-scale error SD will be estimated for each unique
#' reference. These log-scale error SDs will be named following the pattern
#' `sigma_ref_ReferenceID`, where `ReferenceID` is one unique value of
#' `fitdata$Reference`, coerced to character if necessary. This reflects an
#' assumption that all concentration vs. time studies in `fitdata` obey the same
#' underlying PK model, but each study may have a different amount of random
#' measurement error (analogous to meta-regression with fixed effects).
#'
#'  Parameters are estimated using [optimx::optimx()] using the 'L-BFGS-B'
#'  algorithm by default, to enforce lower and upper bounds on parameters.
#'  Parameter starting values are set in [get_starts()]. Parameter lower bounds
#'  are set in [get_lower_bounds()]. Parameter upper bounds are set in
#'  [get_upper_bounds()].
#'
#'  The objective function to be maximized is the log-likelihood, defined in
#'  [log_likelihood()]. It models the residuals as log-normal. The
#'  objective-function gradient is defined analytically in
#'  [grad_log_likelihood()].
#'
#'  Parameters are optimized on the log scale, to reduce scaling issues. This
#'  means that parameter starting values and lower/upper bounds are
#'  log-transformed before being passed to [optimx::optimx()], so if any
#'  starting values or bounds are zero or negative, they will be non-finite when
#'  log-transformed, and [optimx::optimx()] will fail with an error.
#'
#' # Expected variables in \code{fitdata}
#'
#' \describe{\item{\code{Time}}{Time for each data point}
#' \item{\code{Value}}{Concentration for each data point}
#' \item{\code{Dose}}{Dose for each data point} \item{\code{DTXSID}}{DSSTox
#' Substance ID identifying the substance being dosed. This function expects
#' only one unique value in this column; if there are multiple values, it will
#' stop with an error.} \item{\code{Species}}{String identifying the species
#' being dosed. This function expects only one unique value in this column; if
#' there are multiple values, it will stop with an error.}
#' \item{\code{Reference}}{String identifying the reference or study for each
#' data point. Each reference is assumed to have its own residual error standard
#' deviation.} }
#'
#' @param fitdata A \code{data.frame} containing the set of concentration vs.
#'   time data to be fitted. See Details for expected variables. See also
#'   `preprocess_data()`.
#' @param model Which general model should be fit for each chemical. Presently,
#'   only "flat", "1compartment", and "2compartment" are implemented.
#' @param modelfun Either "analytic" or "full" -- whether to fit using the
#'   analytic solution to the model, or the full ODE model. Presently,
#'   "analytic" is recommended (because the analytic solution is exact and much
#'   faster). For `modelfun = 'full'`, only `model = '1compartment'` is
#'   supported.
#' @param pool_sigma Logical: Whether to pool all data (estimate only one error
#'   standard deviation) or not (estimate separate error standard deviations for
#'   each reference). Default FALSE to estimate separate error SDs for each
#'   reference. (If `fitdata` only includes one reference, `pool_sigma` will have
#'   no effect, because only one error SD would be estimated in the first place.)
#' @param LOQ_factor `fitdata$LOQ` will be multiplied by `LOQ_factor` to get the
#'   effective LOQ. Default `LOQ_factor` is 2.
#' @param get_starts_args Any additional arguments to [get_starts()] (other than
#'   `model` and `fitdata`, which are always passed). Default NULL to accept the
#'   default arguments for [get_starts()].
#' @param get_lower_args Any additional arguments to [get_lower_bounds()] (other
#'   than `model` and `fitdata`, which are always passed). Default NULL to
#'   accept the default arguments for [get_lower_bounds()].
#' @param get_upper_args Any additional arguments to [get_upper_bounds()] (other
#'   than `model` and `fitdata`, which are always passed). Default NULL to
#'   accept the default arguments for [get_upper_bounds()].
#' @param optimx_args A named list of additional arguments to
#'   [optimx::optimx()], other than `par`, `fn`, `gr`, `lower`, and `upper`. Default is:
#'
#'    ```
#'     list(
#'           "method" = "L-BFGS-B",
#'           "control" = list("factr" = 1e7,
#'                            "maximize" = TRUE)
#'          )
#'    ```
#'  Note that `method` should allow lower and upper bounds (box constraints),
#'  since these will be supplied. For example, `method` could be "bobyqa" (see
#'  [minqa::bobyqa()]).
#'
#'

analyze_subset <- function(fitdata,
                           model,
                           modelfun,
                           pool_sigma = FALSE,
                           LOQ_factor = 2,
                           get_opt_params_args = NULL,
                           get_starts_args = NULL,
                           get_lower_args = NULL,
                           get_upper_args = NULL,
                           optimx_args = list(
                             "method" = "L-BFGS-B",
                             "control" = list("factr" = 1e7,
                                              "maximize" = TRUE)
                           ),
                           suppress.messages = FALSE,
                           sig.figs = 5,
                           factr = 1e7){

  #convert back to data.frame
  fitdata <- as.data.frame(fitdata)

  this.dtxsid <- unique(fitdata$DTXSID)
  if(length(this.dtxsid) > 1) stop("analyze_subset(): More than one DTXSID in data")

  this.species <- unique(fitdata$Species)
  if(length(this.species) > 1) stop("analyze_subset(): More than one species in data")

  if(!suppress.messages){
    message(paste0("Beginning analysis for:",
                  "Chemical = ",
                   this.dtxsid,
                   "\n",
                   "Species = ",
                   this.species,
                   "\n",
                   "Number of references = ",
                   length(unique(fitdata$Reference)),
                   "\n",
                   "Number of observations = ",
                   nrow(fitdata)
    )
    )
  }

  #get parameter names and
  #determine whether to optimize each of these parameters or not
 par_DF <- do.call(get_opt_params,
                   list("model" = model,
                              "fitdata" = fitdata,
                          "pool_sigma" = pool_sigma,
                        "suppress.messages" = suppress.messages))
#get lower bounds
 par_DF <- do.call(get_lower_bounds,
                   list("par_DF" = par_DF,
                          "model" = model,
                          "fitdata" = fitdata,
                        "pool_sigma" = pool_sigma,
                          "suppress.messages" = suppress.messages))


 #get upper bounds
 par_DF <- do.call(get_upper_bounds,
                   list("par_DF" = par_DF,
                       "model" = model,
                       "fitdata" = fitdata,
                       "pool_sigma" = pool_sigma,
                       "suppress.messages" = suppress.messages))

 #get starting values
 par_DF <- do.call(get_starts,
                   list("par_DF" = par_DF,
                       "model" = model,
                       "fitdata" = fitdata,
                       "pool_sigma" = pool_sigma,
                       "suppress.messages" = suppress.messages))

  #Types of fitted param values to return
  fitted_types <- c("Fitted arithmetic mean",
                    "Fitted arithmetic std dev",
                    "Fitted geometric mean",
                    "Fitted geometric std dev",
                    "Fitted mode",
                    "Fitted log-scale mean",
                    "Fitted log-scale std dev")

  #If the number of parameters is >= the number of data points,
  #then throw back everything NA with a message,
  #because there is no point wasting time trying to fit them.
  if (sum(par_DF$optimize_param) >= nrow(fitdata)){
    #bind together long-form parameter data.frames with a new column "param.value.type"
    #include all of the fitted param.value.types that *would* be provided if we did a fit:
   out_DF <- par_DF
   out_DF[, fitted_types] <- NA_real_

   #check whether there is more than one reference or not
   nref <- length(unique(fitdata$Reference))

    if (nref>1) {
      out_DF$Reference <- paste(sort(unique(fitdata$Reference)),
                                collapse =", ")
      out_DF$Data.Analyzed <- "Joint Analysis"
    } else {
      out_DF$Reference <- unique(fitdata$Reference)
      out_DF$Data.Analyzed <- unique(fitdata$Reference)
    }

    #fill in the loglike and AIC with NA s since no fit was done
    out_DF$LogLikelihood <-  NA_real_
    out_DF$AIC <-  NA_real_

    #include a message about why no fit was done
    msg <- paste("For chemical ", this.dtxsid, " there were ", length(opt_params),
                 " parameters and only ",
                 nrow(fitdata), " data points. Optimization aborted.", sep = "")
    out_DF$message <- msg
    if(!suppress.messages){
    message(msg)
    }
    return(out_DF)
  } #end  if (sum(par_DF$opt.par) >= nrow(fitdata))

  #otherwise, if we have enough data, proceed with the fit.

  #From par_DF, get vectors of:
  #(note all of these are log transformed!!)
  #(parameters are optimized on the log scale,
  #to try and reduce scaling issues)

  #params to be optimized
  opt_params <- log(par_DF[par_DF$optimize_param %in% TRUE,
                       "start_value"])
  names(opt_params) <- par_DF[par_DF$optimize_param %in% TRUE,
                              "param_name"]

  #params to be held constant
  const_params <- log(par_DF[!(par_DF$optimize_param %in% TRUE),
                        "start_value"])
  names(const_params) <- par_DF[!(par_DF$optimize_param %in% TRUE),
                                "param_name"]

  #param lower bounds (only on params to be optimized)
  lower_params <- log(par_DF[par_DF$optimize_param %in% TRUE,
                       "lower_bound"])
  names(lower_params) <- par_DF[par_DF$optimize_param %in% TRUE,
                                "param_name"]

  #param upper bounds (only on params to be optimized)
  upper_params <- log(par_DF[par_DF$optimize_param %in% TRUE,
                         "upper_bound"])
  names(upper_params) <- par_DF[par_DF$optimize_param %in% TRUE,
                                "param_name"]



  all_data_fit <- tryCatch({
    tmp <- do.call(optimx::optimx,
                   args = c(
                     #
                     list(par = opt_params,
                                 fn = log_likelihood,
                          gr = grad_log_likelihood,
                                 #lower = lower_params,
                                 upper = upper_params),
                            #fn, gr, method, and control
                            optimx_args,
                            #... additional args to log_likelihood and grad_log_likelihood
                            list(
                                 const_params = const_params,
                                 DF = fitdata,
                                 modelfun = modelfun,
                                 model = model,
                                 LOQ_factor = LOQ_factor,
                                 force_finite = TRUE)
                            )
                   )
    tmp <- optimx::optimx()
    #tmp is a 1-row data.frame with one variable for each fitted param,
    #plus variables with info on fitting (number of evals, convergence code, etc.)
    #collect any messages from optimx -- in attribute "details" (another data.frame)
    tmp$message <- attr(tmp, "details")[, "message"]
    tmp
  },
           error = function(err){
             #just get the error message
             err$message
           })

  #If fit failed, then return everything as NA and record the message
  if(!is.data.frame(all_data_fit)){
  out_DF <- par_DF
  out_DF[, fitted_types] <- NA_real_

  nref <- length(unique(fitdata$Reference))

  if (nref>1) {
    out_DF$Reference <- paste(sort(unique(fitdata$Reference)),
                              collapse =", ")
    out_DF$Data.Analyzed <- "Joint Analysis"
  } else {
    out_DF$Reference <- unique(fitdata$Reference)
    out_DF$Data.Analyzed <- unique(fitdata$Reference)
  }

  #fill in the loglike and AIC with NA s since no fit was done
  out_DF$LogLikelihood <-  NA_real_
  out_DF$AIC <-  NA_real_

  #include a message about why no fit was done
  msg <- paste("Optimization failed. Error message:",
               all_data_fit)
  out_DF$message <- msg

  if(!suppress.messages){
    message(msg)
  }
  return(out_DF)
  }


#post-processing of fit results

  #Get MLE params
  ln_means <- as.vector(stats::coef(all_data_fit))
  names(ln_means) <- names(opt_params)

  if (!suppress.messages){
    message(paste("Optimized values:  ",
                  paste(apply(data.frame(Names = names(ln_means),
                                         Values = sapply(ln_means,exp),
                                         stringsAsFactors = FALSE),
                              1, function(x) paste(x, collapse = ": ")),
                        collapse = ", "),"\n", sep = ""))

  }

  #Get SDs from Hessian
  #Calculate Hessian using function from numDeriv
  numhess <- numDeriv::hessian(func = objfun,
                               x = ln_means,
                               method = 'Richardson')
  #try inverting Hessian to get SDs
  ln_sds <- tryCatch(diag(solve(numhess)) ^ (1/2),
                     error = function(err){
                       #if hessian can't be inverted
                       if (!suppress.messages) message("Hessian can't be inverted, using pseudovariance matrix to estimate parameter uncertainty.\n")
                       return(diag(chol(MASS::ginv(numhess),
                                        pivot = TRUE)) ^ (1/2)) #pseduovariance matrix
                       #see http://gking.harvard.edu/files/help.pdf
                     })
  names(ln_sds) <- names(ln_means)

  #If any of the SDs are NaN,
  #repeat optimization with smaller convergence tolerance
  while(any(is.nan(ln_sds)) & factr > 1) {

    if (!suppress.messages) message("One or more parameters has NaN standard deviation, repeating optimization with smaller convergence tolerance.\n")
    #then redo optimization
    #bump up starting values for log-scale sigmas by 1
    opt_params[grepl(x = names(opt_params),
                     pattern = "sigma")] <- opt_params[grepl(x = names(opt_params),
                                                              pattern = "sigma")]+1
    #if anything snapped to the bounds,
    #tighten log-scale bounds by 0.1
    opt_params[opt_params <= lower_params] <- lower_params[opt_params <= lower_params] + 0.1
    opt_params[opt_params >= upper_params] <- upper_params[opt_params >= upper_params] - 0.1

    if (!suppress.messages) cat(paste("Initial values:    ", paste(apply(data.frame(Names = names(opt_params),
                                                                                    Values = unlist(lapply(opt_params,exp)),
                                                                                    stringsAsFactors = F),
                                                                         1, function(x) paste(x, collapse = ": ")),
                                                                   collapse = ", "), "\n", sep = ""))
    factr <- factr / 10 #reducing factr by 10 (requiring closer convergence)

    #use general-purpose optimizer to optimize 1-compartment params to fit data
    #optimize by maximizing log-likelihood
    all_data_fit <- optimx::optimx(par = opt_params,
                                   fn = objfun,
                                   # lower=lower,
                                   upper=upper_params,
                                   method = "L-BFGS-B",
                                   hessian = FALSE,
                                   control = list(factr = factr,
                                                  maximize = TRUE))

    ln_means <- as.vector(stats::coef(all_data_fit))
    names(ln_means) <- names(opt_params)
    if (!suppress.messages) message(paste("Optimized values:  ",paste(apply(data.frame(Names = names(ln_means),
                                                                                   Values = sapply(ln_means, exp),
                                                                                   stringsAsFactors=F),
                                                                        1, function(x) paste(x, collapse=": ")),
                                                                  collapse = ", "), "\n", sep = ""))

    #Get SDs from Hessian
    #Calculate Hessian using function from numDeriv
    numhess <- numDeriv::hessian(func = objfun,
                                 x = ln_means,
                                 method = 'Richardson')
    ln_sds <- tryCatch(diag(solve(numhess)) ^ (1/2),
                       error = function(err){
                         #if hessian can't be inverted
                         if (!suppress.messages) cat("Hessian can't be inverted, using pseudovariance matrix to estimate parameter uncertainty.\n")
                         return(diag(chol(MASS::ginv(numhess),
                                          pivot = TRUE)) ^ (1/2)) #pseduovariance matrix
                         #see http://gking.harvard.edu/files/help.pdf
                       })
    names(ln_sds) <- names(ln_means)
  }

  #Arithmetic means:
  #assume log-scale means are for a log-normal distribution
  arith_means <- exp(ln_means + ln_sds ^ 2 / 2)
  names(arith_means) <- names(ln_means)

  #Arithmetic SD
  arith_sd <- (exp(ln_sds ^ 2) - 1)*sqrt(exp(2 * ln_means+
                                             ln_sds ^ 2))

  names(arith_sd) <- names(ln_means)

  #Geometric means
  geo_means <- exp(ln_means)
  names(geo_means) <- names(ln_means)

  #Geometric SDs
  geo_sd <- exp(ln_sds)
  names(geo_sd) <- names(ln_means)

  #Modes
  modes <- exp(ln_means - ln_sds ^ 2)
  names(modes) <- names(ln_means)

  #Produce a data frame of fitted parameters
  fit_DF <- data.frame(arith_means,
                   arith_sd,
                   geo_means,
                   geo_sd,
                   modes,
                   ln_means,
                   ln_sds)
  names(fit_DF) <-  fitted_types
  fit_DF$param_name <- names(ln_means)

  #Merge it with the original data frame of parameters
  out_DF <- merge(par_DF,
                  fit_DF,
                  by = "param_name",
                  all = TRUE)

  nref <- length(unique(fitdata$Reference))

  if (nref>1) {
    out_DF$Reference <- paste(sort(unique(fitdata$Reference)),
                              collapse =", ")
    out_DF$Data.Analyzed <- "Joint Analysis"
  } else {
    out_DF$Reference <- unique(fitdata$Reference)
    out_DF$Data.Analyzed <- unique(fitdata$Reference)
  }

  #Add log-likelihood and AIC values
  out_DF$LogLikelihood <- -all_data_fit$value
  out_DF$AIC <- 2 * length(opt_params) + 2 * all_data_fit$value

  #Check for red flags
  #Initialize flag to blank string...
  out_DF$flag <- ""

  #if anything has not moved from its starting value:
  out_DF[out_DF$optimize_param %in% TRUE &
           out_DF$`Fitted log-scale mean` %in% log(out_DF$start_value),
         "flag"] <- paste0(out_DF[out_DF$optimize_param %in% TRUE &
                                    out_DF$`Fitted log-scale mean` %in%
                                    log(out_DF$start_value),
                                  "flag"],
                           "Fitted log-scale mean equal to log(start value). ")

  #if log-scale std dev is 0 for any parameters
  out_DF[out_DF$optimize_param %in% TRUE &
           out_DF$`Fitted log-scale std dev`%in% 0,
         "flag"] <- paste0(out_DF[out_DF$optimize_param %in% TRUE &
                                    out_DF$`Fitted log-scale std dev`%in% 0,
                                  "flag"],
                           "Fitted log-scale std dev = 0. ")

  #if log-scale std dev is NaN for any parameters
  out_DF[out_DF$optimize_param %in% TRUE &
           !is.finite(out_DF$`Fitted log-scale std dev`),
         "flag"] <- paste0(out_DF[out_DF$optimize_param %in% TRUE &
                                    !is.finite(out_DF$`Fitted log-scale std dev`),
                                  "flag"],
                           "Fitted log-scale std dev is not finite. ")

  #if any parameters are exactly at their bounds
  out_DF[out_DF$optimize_param %in% TRUE &
           (out_DF$`Fitted log-scale mean` %in% log(out_DF$upper_bound) |
           out_DF$`Fitted log-scale mean` %in% log(out_DF$lower_bound)),
         "flag"] <- paste0(out_DF[out_DF$optimize_param %in% TRUE &
                                    (out_DF$`Fitted log-scale mean` %in%
                                       log(out_DF$upper_bound) |
                                       out_DF$`Fitted log-scale mean` %in%
                                       log(out_DF$lower_bound)),
                                  "flag"],
                           "Fitted log-scale mean is equal to lower or upper bound.")

return(out_DF)


}

#' Fit a PK model to one set of concentration vs. time data
#'
#' Fit a specified PK model to one set of concentration vs. time data and return
#' a set of fitted parameter values.
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
#' # Specifications for `fun_start`
#'
#' ## Arguments
#'
#' `fun_start` must accept the same named arguments as \code{get_starts},
#' i.e. the parameters \code{fitdata}, \code{model}. It may
#' accept additional arguments, to be specified in a named list in
#' \code{fun_start_args}.
#'
#' ## Return value
#'
#' `fun_start` must return a \code{data.frame} of model parameters with the
#' following variables: \describe{ \item{\code{param.name}}{Character: the names
#' of the PK model parameters. Must include all of the parameter names expected
#' by the model function defined by \code{model}, as well as
#' all of the reference-specific error standard deviations, named as
#' \code{paste("sigma", unique(fitdata$Reference), sep = ".)}}.
#' \item{\code{param.value}}{Numeric: A starting value for each of the PK model
#' parameters.} \item{\code{param.value.source}{Character: A message describing
#' briefly how the starting value was calculated for each parameter. This may be
#' blank or NA if you do not wish to supply any message, but the variable must
#' be present.}} }
#'
#' ## Example
#'
#' If you want to accept the default starting values for most parameters, but
#' you want to change the starting values for the reference-specific error SDs
#' so that they are all fixed at some constant value (e.g., 1), you could define
#' your own function as follows:
#'
#' ```{r, eval = FALSE}
#' my_starts <- function(fitdata,
#'                        model,
#'                        sigma_start_const = 1){
#'    #Get default starting values for all params
#'    starts_df <- get_starts(fitdata, model)
#'    #Overwrite defaults for reference-specific error SDs
#'    starts_df[grepl(x = starts_df$param.name,
#'                    pattern = "sigma"),
#'                    c("param.value",
#'                    "param.value.source")] <- list(sigma_start_const,
#'                    "Set to constant")
#'    return(starts_df)
#' }
#' ```
#'
#' Then when calling `analyze_subset`, you would specify `fun_start = my_starts`
#' and `fun_start_args = list(sigma_start_const=1)`.
#'
#' # Specifications for `fun_lower`
#'
#' ## Arguments
#'
#' `fun_lower` must accept the same named arguments as
#' \code{get_lower_bounds}, i.e. the parameters `par_DF`, `fitdata`, `model`. It
#' may accept additional arguments, to be specified in a named list in
#' \code{fun_lower_args}.
#'
#' ## Return value
#'
#' `fun_lower` must return a \code{data.frame} of model parameters with the
#' following variables: \describe{ \item{\code{param.name}}{Character: the names
#' of the PK model parameters. Must include all of the parameter names expected
#' by the model function defined by \code{model}, as well as
#' all of the reference-specific error standard deviations, named as
#' \code{paste("sigma", unique(fitdata$Reference), sep = ".)}}.
#' \item{\code{param.lower}}{Numeric: A lower bound for each of the PK model
#' parameters.} \item{\code{param.lower.source}{Character: A message describing
#' briefly how the lower-bound value was calculated for each parameter. This may
#' be blank or NA if you do not wish to supply any message, but the variable
#' must be present.}} }
#'
#' # Specifications for `fun_upper`
#'
#' ## Arguments
#'
#' `fun_upper` must accept the same named arguments as
#' \code{get_upper_bounds}, i.e. the parameters \code{fitdata}, \code{model}. It
#' may accept additional arguments, to be specified in a named list in
#' \code{fun_upper_args}.
#'
#' ## Return value
#'
#' `fun_upper` must return a \code{data.frame} of model parameters with the
#' following variables: \describe{ \item{\code{param.name}}{Character: the names
#' of the PK model parameters. Must include all of the parameter names expected
#' by the model function defined by \code{model}, as well as
#' all of the reference-specific error standard deviations, named as
#' \code{paste("sigma", unique(fitdata$Reference), sep = ".)}}.
#' \item{\code{param.upper}}{Numeric: An upper bound for each of the PK model
#' parameters.} \item{\code{param.upper.source}{Character: A message describing
#' briefly how the upper-bound value was calculated for each parameter. This may
#' be blank or NA if you do not wish to supply any message, but the variable
#' must be present.}} }
#'
#' @param fitdata A \code{data.frame} containing the set of concentration vs.
#'   time data to be fitted. See Details for expected variables. See also
#'   \code{\link{preprocess_data}}.
#' @param model Which general model should be fit for each chemical. Presently,
#'   only "1compartment" and "2compartment" are implemented.
#' @param modelfun Either "analytic" or "full" -- whether to fit using the
#'   analytic solution to the model, or the full ODE model. Presently,
#'   "analytic" is recommended (because the analytic solution is exact and much
#'   faster).
#' @param fun_start Function to calculate starting values for each parameter to
#'   be fitted. Default \code{get_starts}. See Details for expected arguments
#'   and return value for this function.
#' @param fun_lower Function to calculate lower bounds for each parameter to be
#'   fitted. Default \code{get_lower_bounds}. See Details for expected arguments
#'   and return value for this function.
#' @param fun_upper Function to calculate upper bounds for each parameter to be
#'   fitted. Default \code{get_upper_bounds}. See Details for expected arguments
#'   and return value for this function.
#' @param fun_start_args A named list of additional arguments to
#'   `fun_start`, if any. Default NULL.
#' @param fun_lower_args A named list of additional arguments to
#'   `fun_lower`, if any. Default NULL.
#' @param fun_upper_args A named list of additional arguments to
#'   `fun_upper`, if any. Default NULL.

analyze_subset <- function(fitdata,
                           model,
                           modelfun,
                           sigma_ref = TRUE,
                           fun_start = get_starts,
                           fun_lower = get_lower_bounds,
                           fun_upper = get_upper_bounds,
                           fun_start_args = NULL,
                           fun_lower_args = NULL,
                           fun_upper_args = NULL,
                           this.reference = NULL,
                           suppress.messages = FALSE,
                           sig.figs = 5){

  this.dtxsid <- unique(fitdata$DTXSID)
  if(length(this.dtxsid) > 1) stop("analyze_subset(): More than one DTXSID in data")

  this.species <- unique(fitdata$species)
  if(length(this.species) > 1) stop("analyze_subset(): More than one species in data")

  #get parameter names and
  #determine whether to optimize each of these parameters or not
 par_DF <- get_opt_params(model = model,
                              fitdata = fitdata,
                          sigma_ref = sigma_ref)
#get lower bounds
 par_DF <- do.call(fun_lower,
                   c(list("par_DF" = par_DF,
                          "model" = model,
                          "fitdata" = fitdata),
                     fun_lower_args))

 #get upper bounds
 par_DF <- do.call(fun_upper,
                   c(list("par_DF" = par_DF,
                          "model" = model,
                          "fitdata" = fitdata),
                     fun_upper_args))

 #get starting values
 par_DF <- do.call(fun_start,
                   c(list("par_DF" = par_DF,
                          "model" = model,
                          "fitdata" = fitdata),
                     fun_start_args))

  #Types of fitted param values to return
  fitted_types <- c("Fitted arithmetic mean",
                    "Fitted arithmetic std dev",
                    "Fitted geometric mean",
                    "Fitted geometric std dev",
                    "Fitted mode")

  #


  #If the number of parameters is >= the number of data points,
  #then throw back everything NA with a message,
  #because there is no point wasting time trying to fit them.
  if (sum(par_DF$optimize_param) >= nrow(fitdata)){
    #bind together long-form parameter data.frames with a new column "param.value.type"
    #include all of the fitted param.value.types that *would* be provided if we did a fit:
   out_DF <- par_DF
   out_DF[, fitted_types] <- NA_real_

    if (is.null(this.reference)) {
      out_DF$Reference <- paste(sort(unique(fitdata$Reference)),
                                collapse =", ")
      out_DF$Data.Analyzed <- Reference
      out_DF[regexpr(",", out_DF$Data.Analyzed) != -1, "Data.Analyzed"] <-  'Joint Analysis'
    } else {
      out_DF$Data.Analyzed <- this.reference
    }

    #fill in the loglike and AIC with NA s since no fit was done
    out_DF$LogLikelihood <-  NA_real_
    out_DF$AIC <-  NA_real_

    #include a message about why no fit was done
    msg <- paste("For chemical ", this.dtxsid, " there were ", length(opt_params),
                 " parameters and only ",
                 nrow(fitdata), " data points. Optimization aborted.", sep = "")
    out_DF$message <- msg

    message(msg)
    return(out_DF)
  } #end  if (sum(par_DF$opt.par) >= nrow(fitdata))

  #otherwise, if we have enough data, proceed with the fit.

  #From par_DF, get vectors of:
  #(note all of these are log transformed!!)

  #params to be optimized
  opt_params <- log(par_DF[par_DF$optimize_param %in% TRUE,
                       "start_value"])
  names(opt_params) <- par_DF$param_name

  #params to be held constant
  const_params <- log(par_DF[!(par_DF$optimize_param %in% TRUE),
                        "start_value"])
  names(const_params) <- par_DF$param_name

  #param lower bounds (only on params to be optimized)
  lower_params <- log(par_DF[par_DF$optimize_param %in% TRUE,
                       "lower_value"])
  names(lower_params) <- par_DF$param_name

  #param upper bounds (only on params to be optimized)
  upper_params <- log(par_DF[par_DF$optimize_param %in% TRUE,
                         "upper_value"])
  names(upper_params) <- par_DF$param_name


#objective function
  objfun <- function(x) {
    foo <- -log_likelihood(params = c(x,
                                      const_params),
                           DT = fitdata,
                           modelfun = modelfun,
                           model = model)
    #method L-BFGS-B requires finite values be returned,
    #so if log-likelihood is NA or -Inf,
    #just return a large negative log.likelihood
    #(= a large positive -log.likelihood)
    if (!is.finite(foo)) foo <- 99999
    return(foo)
  }



  all_data_fit <- tryCatch({
    tmp <- optimx::optimx(par = log(opt_params),
                                          fn = objfun,
                                          #lower = lower_params,
                                          upper = upper_params,
                                          method = "L-BFGS-B", #box constraints
                                          hessian = FALSE,
                                          control = list(factr = factr))
    #collect any messages from optimx
    tmp$message <- attr(tmp, "details")$message
  },
           error = function(err){
             data.frame(message = err$message)

           })

#post-processing of fit results

  #Get MLE params
  ln_means <- as.vector(stats::coef(all_data_fit))
  names(ln_means) <- names(opt_params)

  if (!suppress.messages) message(paste("Optimized values:  ", paste(apply(data.frame(Names = names(ln_means),
                                                                                  Values = sapply(ln_means,exp),
                                                                                  stringsAsFactors = FALSE),
                                                                       1, function(x) paste(x, collapse = ": ")),
                                                                 collapse = ", "),"\n", sep = ""))

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
    #bump up log-scale sigmas by 1
    #(equivalent to multiplying sigma by 10)
    opt_params[regexpr("sigma", names(opt_params)) != -1] <- opt_params[regexpr("sigma", names(opt_params))!=-1]+1
    #if anything snapped to the bounds,
    #tighten log-scale bounds by 0.1
    #(equivalent to tightening by a factor of 10)
    opt_params[opt_params <= lower] <- lower[opt_params <= lower] + 0.1
    opt_params[opt_params >= upper] <- upper[opt_params >= upper] - 0.1

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
                                   upper=upper,
                                   method = "L-BFGS-B",
                                   hessian = FALSE,
                                   control = list(factr = factr))

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
                   modes)
  names(fit_DF) <-  fitted_types
  fit_DF[["Fitted log-scale mean"]] <- ln_means
  fit_DF[["Fitted log-scale std dev"]] <- ln_sds
  fit_DF$param_name <- names(ln_means)

  #Merge it with the original data frame of parameters
  out_DF <- merge(par_DF,
                  fit_DF,
                  by = "param_name",
                  all = TRUE)
  #Add log-likelihood and AIC values
  out_DF$LogLikelihood <- -all_data_fit$value
  out_DF$AIC <- 2 * length(opt_params) + 2 * all_data_fit$value

  #Check for red flags
  #if anything has not moved from its starting value:
  out_DF[out_DF$optimize_param %in% TRUE &
           out_DF$`Fitted log-scale mean` == log(out_DF$start_value),
         "flag"] <- "Fitted log-scale mean equal to log(start value)"

  #if log-scale std dev is 0 for any parameters
  out_DF[out_DF$optimize_param %in% TRUE &
           out_DF$`Fitted log-scale std dev`==0,
         "flag"] <- "Fitted log-scale std dev = 0"

  #if log-scale std dev is NaN for any parameters
  out_DF[out_DF$optimize_param %in% TRUE &
           !is.finite(out_DF$`Fitted log-scale std dev`),
         "flag"] <- "Fitted log-scale std dev is not finite"

  #if any parameters are exactly at their bounds
  out_DF[out_DF$optimize_param %in% TRUE &
           (out_DF$`Fitted log-scale mean` == log(out_DF$upper_value) |
           out_DF$`Fitted log-scale mean` == log(out_DF$lower_value)),
         "flag"] <- "Fitted log-scale mean is equal to lower or upper bound"


return(out_DF)


}

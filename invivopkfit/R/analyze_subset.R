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
#'  `kelim`, and possibly `kgutabs` and `Fgutabs` or `Fgutabs_Vdist` (see below).
#'  - for `model = '2compartment'`, the parameters estimated are `V1`, `kelim`,
#'  `k12`, `k21`, and possibly `kgutabs` and `Fgutabs` or `Fgutabs_Vdist`.
#'
#'  For 1-compartment and 2-compartment models, `kgutabs` will be estimated from
#'  the data only if `fitdata` includes oral dosing data; otherwise it will be
#'  set to NA.
#'
#'  `Fgutabs` will be estimated from the data only if `fitdata` includes both
#'  oral and IV data. If `fitdata` includes oral data but not IV data, then
#'  `Fgutabs_Vdist` (1-compartment) or `Fgutabs_V1` (2-compartment) will be
#'  estimated instead of `Fgutabs` and `Vdist` or `V1` (because only the ratio
#'  of `Fgutabs` and `Vdist` is identifiable in that case).
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
#'  # Parameter estimation via numerical optimization
#'
#'  Parameters are estimated using [optimx::optimx()] by default using the
#'  [minqa::bobyqa()] algorithm. (The algorithm is specified in
#'  `optimx_args$method`.) The objective function to be maximized is the
#'  log-likelihood, defined in [log_likelihood()]. It models the residuals as
#'  log-normal. Note that the log-likelihood is *maximized*, not minimized.
#'  Maximization is specified by setting argument `optimx_args$control$maximize
#'  = TRUE`. Parameter starting values are set in [get_starts()]. Parameter
#'  lower bounds are set in [get_lower_bounds()]. Parameter upper bounds are set
#'  in [get_upper_bounds()].
#'
#'  By default, the objective function gradient is estimated numerically.
#'
#'  ## Parameter standard deviations via numerical Hessian
#'
#'  Parameter standard deviations (uncertainty) are estimated as the square root
#'  of the diagonal of the inverse of the Hessian (the matrix of second
#'  derivatives of the objective function). The inverse Hessian approximates the
#'  variance-covariance matrix of estimated parameters. The square root of its
#'  diagonal approximates parameter standard deviations.
#'
#'  If the numerically-estimated Hessian cannot be inverted, a
#'  pseudoinverse is attempted.
#'
#'  If the square root of the inverse Hessian contains NaN, then if
#'  `optimx_args$method %in% 'L-BFGS-B'`, optimization is repeated with a
#'  smaller convergence tolerance, until either there are no more NaNs or
#'  convergence tolerance factor is equal to 1 (i.e. tolerance is at machine
#'  epsilon). This is not done if `optimx_args$method %in% 'bobyqa'` because
#'  that method does not use a convergence tolerance control parameter.
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
#'   [optimx::optimx()], other than `par`, `fn`, `lower`, and `upper`. Default is:
#'
#'    ```
#'     list(
#'           "method" = "bobyqa",
#'           "itnmax" = 1e6,
#'           "control" = list("maximize" = TRUE)
#'          )
#'    ```
#'  See documentation for [optimx::optimx()] for arguments and details. Note
#'  lower and upper bounds (box constraints) will be supplied; if you want them
#'  to be respected, please choose a method that allows box constraints (e.g.
#'  "bobyqa" or "L-BFGS-B").
#'
#'

analyze_subset <- function(fitdata,
                           model,
                           modelfun,
                           pool_sigma = FALSE,
                           get_starts_args = NULL,
                           get_lower_args = NULL,
                           get_upper_args = NULL,
                           optimx_args = list(
                             "method" = "bobyqa",
                             "itnmax" = 1e6,
                             "control" = list("maximize" = TRUE,
                                              "kkt" = FALSE)
                           ),
                           suppress.messages = FALSE){

  #convert back to data.frame
  fitdata <- as.data.frame(fitdata)

  this.dtxsid <- unique(fitdata$DTXSID)
  if(length(this.dtxsid) > 1) stop("analyze_subset(): More than one DTXSID in data")

  this.species <- unique(fitdata$Species)
  if(length(this.species) > 1) stop("analyze_subset(): More than one species in data")

  if(!suppress.messages){
    message(paste0("Beginning analysis for:\n",
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
                   nrow(fitdata),
                  " (iv: ",
                  sum(fitdata$Route %in% "iv"),
                  "; po: ",
                  sum(fitdata$Route %in% "po"),
                  ")"
    )
    )
  }

  #assign default max iterations/function evals if missing
  if(is.null(optimx_args$itnmax)){
    optimx_args$itnmax <- max(5e4,
                      5000*round(
                        sqrt(
                          length(opt_params)+1
                        )
                      )
        )


    if(!suppress.messages){
      message(paste("optimx_args does not include item 'itnmax'.",
                    "Setting optimx_args$itnmax =",
                    "max(5e4, 5000*round(sqrt(n+1)))",
                    "where n = number of params to be optimized.",
                    "(This is the optimx() default)"))
    }
  }

  #assign convergence tolerance factor if missing and relevant
  if(optimx_args$method %in% "L-BFGS-B" &
     is.null(optimx_args$control$factr)){
    optimx_args$control$factr <- 1e7
    if(!suppress.messages){
      message(paste("Optimization method is",
                    optimx_args$method,
                    "but control argument factr was not provided;",
                    "setting factr = 1e7"))
    }
  }

  #get parameter names and
  #determine whether to optimize each of these parameters or not
 par_DF <- do.call(get_opt_params,
                   list("model" = model,
                              "fitdata" = fitdata,
                          "pool_sigma" = pool_sigma,
                        "suppress.messages" = suppress.messages))

 #Types of fitted param values to return
 fitted_types <- c("Fitted mean",
                   "Fitted std dev")

 #If the number of parameters is >= the number of detected data points,
 #then throw back everything NA with a message,
 #because there is no point wasting time trying to fit them.
 if (sum(par_DF$optimize_param) >= sum(!is.na(fitdata$Value))){

   out_DF <- par_DF[par_DF$optimize_param %in% TRUE, ]

   #add NA for lower, upper bounds and start values
   out_DF[, c("lower_bound",
              "lower_bound_msg",
              "upper_bound",
              "upper_bound_msg",
              "start_value",
              "start_value_msg")] <- list(NA_real_,
                                          NA_character_,
                                          NA_real_,
                                          NA_character_,
                                          NA_real_,
                                          NA_character_)

   out_DF[, fitted_types] <- NA_real_

   #check whether there is more than one reference or not
   nref <- length(unique(fitdata$Reference))

   if (nref>1) {
     if(pool_sigma %in% FALSE){
       out_DF$Reference <- paste(sort(unique(fitdata$Reference)),
                                 collapse =", ")
       out_DF$Data.Analyzed <- "Joint Analysis"
     }else{
       out_DF$Reference <- paste(sort(unique(fitdata$Reference)),
                                 collapse =", ")
       out_DF$Data.Analyzed <- "Pooled Analysis"
     }
   }else {
     out_DF$Reference <- unique(fitdata$Reference)
     out_DF$Data.Analyzed <- "Single-Reference Analysis"
   }

   #fill in the loglike and AIC with NA s since no fit was done
   out_DF$LogLikelihood <-  NA_real_
   out_DF$AIC <-  NA_real_

   #include a message about why no fit was done
   msg <- paste("For chemical ", this.dtxsid, " there were ",
                sum(par_DF$optimize_param),
                " parameters to be estimated (",
                paste(par_DF[par_DF$optimize_param %in% TRUE, "param_name"],
                      collapse = ", "),
                ") and only ",
                sum(!is.na(fitdata$Value)),
                " detected data points. Optimization aborted.",
                sep = "")
   out_DF$message <- msg
   out_DF$flag <- NA_character_
   #Record the unique routes in this dataset
   #Route info provides context for why some parameters were/were not estimated
   out_DF$Routes <- paste("iv: ",
                          sum(fitdata$Route %in% "iv"),
                          "; po: ",
                          sum(fitdata$Route %in% "po"))
   out_DF$fevals <- NA_integer_
   out_DF$convcode <- NA_integer_
   out_DF$niter <- NA_integer_
   out_DF$method <- optimx_args$method
   #record control params
   out_DF[paste0("control_",
                 names(optimx_args$control))] <- optimx_args$control
   if(!suppress.messages){
     message(msg)
   }
   return(out_DF)
 } #end  if (sum(par_DF$opt.par) >= nrow(fitdata))

 #else -- continue
#get lower bounds
 par_DF <- do.call(get_lower_bounds,
                   c(list("par_DF" = par_DF,
                          "model" = model,
                          "fitdata" = fitdata,
                        "pool_sigma" = pool_sigma,
                          "suppress.messages" = suppress.messages),
                     get_lower_args))


 #get upper bounds
 par_DF <- do.call(get_upper_bounds,
                   c(list("par_DF" = par_DF,
                       "model" = model,
                       "fitdata" = fitdata,
                       "pool_sigma" = pool_sigma,
                       "suppress.messages" = suppress.messages),
                     get_upper_args)
 )

 #get starting values
 par_DF <- do.call(get_starts,
                   c(list("par_DF" = par_DF,
                       "model" = model,
                       "fitdata" = fitdata,
                       "pool_sigma" = pool_sigma,
                       "suppress.messages" = suppress.messages),
                     get_starts_args))

 #Initialize out_DF
 #There will be one row for each parameter
out_DF <- par_DF[par_DF$optimize_param %in% TRUE, ]

  #From par_DF, get vectors of:

  #params to be optimized
  opt_params <- par_DF[par_DF$optimize_param %in% TRUE,
                       "start_value"]
  names(opt_params) <- par_DF[par_DF$optimize_param %in% TRUE,
                              "param_name"]

  #param lower bounds (only on params to be optimized)
  lower_params <- par_DF[par_DF$optimize_param %in% TRUE,
                         "lower_bound"]
  names(lower_params) <- par_DF[par_DF$optimize_param %in% TRUE,
                                "param_name"]

  #param upper bounds (only on params to be optimized)
  upper_params <- par_DF[par_DF$optimize_param %in% TRUE,
                         "upper_bound"]
  names(upper_params) <- par_DF[par_DF$optimize_param %in% TRUE,
                                "param_name"]

  #get constant params if any
  if(any(par_DF$optimize_param %in% FALSE &
         par_DF$use_param %in% TRUE)){
  const_params <- par_DF[par_DF$optimize_param %in% FALSE &
                           par_DF$use_param %in% TRUE,
                         "start_value"]

  names(const_params) <- par_DF[par_DF$optimize_param %in% FALSE &
                                  par_DF$use_param %in% TRUE,
                                "param_name"]
  }else{
    const_params <- NULL
  }

  if (!suppress.messages){
    message(paste0("Estimating parameters ",
                   paste(names(opt_params), collapse = ", "),
                   "\n",
    "Using optimx::optimx(), ",
                   "method = ", optimx_args$method, "\n",
                   "Convergence tolerance factr = ",
                   optimx_args$control$factr, "\n",
    "..."))

    message(paste("Initial values:    ",
                  paste(apply(data.frame(Names = names(opt_params),
                                         Values = unlist(opt_params),
                                         stringsAsFactors = F),
                              1, function(x) paste(x, collapse = ": ")),
                        collapse = ", "),
                  sep = ""))
  }

  all_data_fit <- tryCatch({
    if(suppress.messages %in% TRUE){
      junk <- capture.output(
        tmp <- do.call(
          optimx::optimx,
          args = c(
            #
            list(par = opt_params,
                 fn = log_likelihood,
                 lower = lower_params,
                 upper = upper_params),
            #method and control
            optimx_args,
            #... additional args to log_likelihood
            list(
              const_params = const_params,
              DF = fitdata,
              modelfun = modelfun,
              model = model,
              force_finite = FALSE
            ) #end list()
          ) #end args = c()
        ) #end do.call
      ) #end capture.output
    }else{
      tmp <- do.call(
        optimx::optimx,
        args = c(
          #
          list(par = opt_params,
               fn = log_likelihood,
               lower = lower_params,
               upper = upper_params),
          #method, and control
          optimx_args,
          #... additional args to log_likelihood
          list(
            const_params = const_params,
            DF = fitdata,
            modelfun = modelfun,
            model = model,
            force_finite = FALSE
          )
        )
      ) #end do.call
    }
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
  out_DF[, fitted_types] <- NA_real_

  nref <- length(unique(fitdata$Reference))

  if (nref>1) {
    if(pool_sigma %in% FALSE){
    out_DF$Reference <- paste(sort(unique(fitdata$Reference)),
                              collapse =", ")
    out_DF$Data.Analyzed <- "Joint Analysis"
    }else{
      out_DF$Reference <- paste(sort(unique(fitdata$Reference)),
                                collapse =", ")
      out_DF$Data.Analyzed <- "Pooled Analysis"
    }
  } else {
    out_DF$Reference <- unique(fitdata$Reference)
    out_DF$Data.Analyzed <- "Single-Reference Analysis"
  }

  #fill in the loglike and AIC with NA s since no fit was done
  out_DF$LogLikelihood <-  NA_real_
  out_DF$AIC <-  NA_real_

  #include a message about why no fit was done
  msg <- paste("Optimization failed. Error message from optimx():",
               all_data_fit)
  out_DF$message <- msg

  out_DF$flag <- NA_character_
  #Record the unique routes in this dataset
  #Route info provides context for why some parameters were/were not estimated
  out_DF$Routes <- paste("iv: ",
                         sum(fitdata$Route %in% "iv"),
                         "; po: ",
                         sum(fitdata$Route %in% "po"))

  out_DF$fevals <- NA_integer_
  out_DF$convcode <- NA_integer_
  out_DF$niter <- NA_integer_
  out_DF$method <- optimx_args$method
  #record control params
  out_DF[paste0("control_",
                names(optimx_args$control))] <- optimx_args$control

  if(!suppress.messages){
    message(msg)
  }
  return(out_DF)
  }


#post-processing of fit results

  #Get MLE params
  means <- as.vector(stats::coef(all_data_fit))
  names(means) <- names(opt_params)

  if (!suppress.messages){
    message(paste("Optimized values:  ",
                  paste(apply(data.frame(Names = names(means),
                                         Values = means,
                                         stringsAsFactors = FALSE),
                              1, function(x) paste(x, collapse = ": ")),
                        collapse = ", "),"\n", sep = ""))

  }

  #Get SDs from Hessian
  #Calculate Hessian using function from numDeriv
  numhess <- numDeriv::hessian(func = function(x){
    -1 * log_likelihood(x,
                        const_params = const_params,
                        DF = fitdata,
                        modelfun = modelfun,
                        model = model,
                        force_finite = FALSE )
  },
  x = means,
  method = 'Richardson')

  #try inverting Hessian to get SDs
  sds <- tryCatch(diag(solve(numhess)) ^ (1/2),
                     error = function(err){
                       #if hessian can't be inverted
                       if (!suppress.messages) {
                         message(paste0("Hessian can't be inverted, ",
                         "using pseudovariance matrix ",
                         "to estimate parameter uncertainty."))
                       }
                       suppressWarnings(tmp <- diag(chol(MASS::ginv(numhess),
                                        pivot = TRUE)) ^ (1/2))
                       return(tmp) #pseudovariance matrix
                       #see http://gking.harvard.edu/files/help.pdf
                     })
  names(sds) <- names(means)

  #Only if
  #If any of the SDs are NaN,
  #repeat optimization with smaller convergence tolerance
  #-- only if method is "L-BFGS-B" -- other methods do not use this
  if(optimx_args$method %in% "L-BFGS-B"){
  while(any(is.nan(sds)) & optimx_args$control$factr > 1) {


    #then redo optimization

    #start from fitted values

    #
    # #update the starting points in par_DF so they will be recorded
    par_DF[match(names(means),
                 par_DF$param_name),
           "start_value"] <- means

    #get new starting values to use
    opt_params <- par_DF[par_DF$optimize_param %in% TRUE,
                             "start_value"]
    names(opt_params) <- par_DF[par_DF$optimize_param %in% TRUE,
                                "param_name"]


    optimx_args$control$factr <- optimx_args$control$factr / 10 #reducing factr by 10 (requiring closer convergence)

    if (!suppress.messages){
      message(paste0("One or more parameters has NaN standard deviation. ",
                     "Repeating optimization with smaller convergence tolerance, ",
                     "starting from result of previous optimization:\n",
                     paste(
                       apply(
                         data.frame(Names = names(opt_params),
                                    Values = unlist(opt_params),
                                    stringsAsFactors = F),
                         1,
                         function(x) paste(x, collapse = ": ")
                       ),
                       collapse = ", "
                     )
      )
      )
      message(paste0("Estimating parameters ",
                     paste(names(opt_params), collapse = ", "),
                     "\n",
                     "Using optimx::optimx(), ",
                     "method = ", optimx_args$method, "\n",
                     "Convergence tolerance factr = ",
                     optimx_args$control$factr, "\n",
                     "..."))

      all_data_fit <-  do.call(optimx::optimx,
                               args = c(
                                 #
                                 list(par = opt_params,
                                      fn = log_likelihood,
                                      lower = lower_params,
                                      upper = upper_params),
                                 #method and control
                                 optimx_args,
                                 list(
                                   const_params = const_params,
                                   DF = fitdata,
                                   modelfun = modelfun,
                                   model = model,
                                   force_finite = FALSE
                                 )
                               )
      )
    }else{

      junk <- capture.output(
        all_data_fit <-  do.call(
          optimx::optimx,
          args = c(
            #
            list(par = opt_params,
                 fn = log_likelihood,
                 lower = lower_params,
                 upper = upper_params),
            #method and control
            optimx_args,
            list(
              const_params = const_params,
              DF = fitdata,
              modelfun = modelfun,
              model = model,
              force_finite = FALSE
            )
          )
        )
      )
    }

    means <- as.vector(stats::coef(all_data_fit))
    names(means) <- names(opt_params)
    if (!suppress.messages){
      message(paste("Optimized values:  ",
                    paste(apply(data.frame(Names = names(means),
                                           Values = means,
                                           stringsAsFactors=F),
                                1, function(x) paste(x, collapse=": ")),
                          collapse = ", "),
                    "\n", sep = ""))
    }

    #Get SDs from Hessian
    #Calculate Hessian using function from numDeriv
    #but use NEGATIVE log likelihood
    numhess <- numDeriv::hessian(func = function(x){
      -1 * log_likelihood(x,
                          const_params = const_params,
                          DF = fitdata,
                          modelfun = modelfun,
                          model = model,
                          force_finite = FALSE )
    },
    x = means,
    method = 'Richardson')

    sds <- tryCatch(diag(solve(numhess)) ^ (1/2),
                       error = function(err){
                         #if hessian can't be inverted
                         if (!suppress.messages){
                           message("Hessian can't be inverted, using pseudovariance matrix to estimate parameter uncertainty.")
                         }
                         suppressWarnings(tmp <- diag(chol(MASS::ginv(numhess),
                                          pivot = TRUE)) ^ (1/2))
                         return(tmp) #pseduovariance matrix
                         #see http://gking.harvard.edu/files/help.pdf
                       })
    names(sds) <- names(means)
  } #end while(any(is.nan(sds)) & optimx_args$control$factr > 1)
  } #end if method %in% "L-BFGS-B"

  #Produce a data frame of fitted parameters
  fit_DF <- data.frame(means,
                   sds)
  names(fit_DF) <-  fitted_types
  fit_DF$param_name <- names(means)

  #Merge it with the original data frame of parameters
  #keep only the ones that were actually fit
  out_DF <- merge(par_DF,
                  fit_DF,
                  by = "param_name",
                  all.y = TRUE,
                  all.x = FALSE)

  nref <- length(unique(fitdata$Reference))

  if (nref>1) {
    if(pool_sigma %in% FALSE){
      out_DF$Reference <- paste(sort(unique(fitdata$Reference)),
                                collapse =", ")
      out_DF$Data.Analyzed <- "Joint Analysis"
    }else{
      out_DF$Reference <- paste(sort(unique(fitdata$Reference)),
                                collapse =", ")
      out_DF$Data.Analyzed <- "Pooled Analysis"
    }
  } else {
    out_DF$Reference <- unique(fitdata$Reference)
    out_DF$Data.Analyzed <- "Single-Reference Analysis"
  }

  #Add log-likelihood and AIC values
  out_DF$LogLikelihood <- all_data_fit$value
  out_DF$AIC <- 2 * length(opt_params) - 2 * all_data_fit$value

  #Check for red flags
  #Initialize flag to blank string...
  out_DF$flag <- ""

  #if anything has not moved from its starting value:
  out_DF[out_DF$optimize_param %in% TRUE &
           out_DF$`Fitted mean` %in% out_DF$start_value,
         "flag"] <- paste0(out_DF[out_DF$optimize_param %in% TRUE &
                                    out_DF$`Fitted mean` %in%
                                    log(out_DF$start_value),
                                  "flag"],
                           "Fitted mean equal to start value. ")

  #if log-scale std dev is 0 for any parameters
  out_DF[out_DF$optimize_param %in% TRUE &
           out_DF$`Fitted std dev`%in% 0,
         "flag"] <- paste0(out_DF[out_DF$optimize_param %in% TRUE &
                                    out_DF$`Fitted std dev`%in% 0,
                                  "flag"],
                           "Fitted std dev = 0. ")

  #if log-scale std dev is NaN for any parameters
  out_DF[out_DF$optimize_param %in% TRUE &
           !is.finite(out_DF$`Fitted std dev`),
         "flag"] <- paste0(out_DF[out_DF$optimize_param %in% TRUE &
                                    !is.finite(out_DF$`Fitted std dev`),
                                  "flag"],
                           "Fitted std dev is not finite. ")

  #if any parameters are exactly at their bounds
  out_DF[out_DF$optimize_param %in% TRUE &
           (out_DF$`Fitted mean` %in% out_DF$upper_bound |
           out_DF$`Fitted mean` %in% out_DF$lower_bound),
         "flag"] <- paste0(out_DF[out_DF$optimize_param %in% TRUE &
                                    (out_DF$`Fitted mean` %in%
                                       out_DF$upper_bound |
                                       out_DF$`Fitted mean` %in%
                                       out_DF$lower_bound),
                                  "flag"],
                           "Fitted mean is equal to lower or upper bound.")

  out_DF[!nzchar(out_DF$flag), "flag"] <- NA_character_


  #Record the unique routes in this dataset
  #Route info provides context for why some parameters were/were not estimated
  out_DF$Routes <- paste("iv: ",
                         sum(fitdata$Route %in% "iv"),
                         "; po: ",
                         sum(fitdata$Route %in% "po"))
  out_DF$message <- "Optimization successful."


  out_DF$fevals <- as.integer(all_data_fit$fevals)
  out_DF$convcode <- as.integer(all_data_fit$convcode)
  out_DF$niter <- as.integer(all_data_fit$niter)
  out_DF$method <- optimx_args$method
  #record control params
  out_DF[paste0("control_",
                names(optimx_args$control))] <- optimx_args$control


return(out_DF)


}

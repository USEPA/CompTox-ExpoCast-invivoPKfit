#'Fit a PK model to one set of concentration vs. time data
#'
#'Fit a specified PK model to one set of concentration vs. time data and return
#'a set of fitted parameter values.
#'
#'Typically this function is not called directly by the user, but is called from
#'[fit_all()], the main fitting function.
#'
#'This function estimates one set of parameters for a specified PK model, based
#'on a set of concentration vs. time data in `fitdata`. `fitdata` is typically
#'data for one chemical and one species, but it may include data from more than
#'one study and more than one dosing route.
#'
#'The parameters to be estimated are those defined in [get_opt_params()] for the
#'model specified in `model` and `modelfun`:
#'
#'  - for `model = 'flat'`, the parameter estimated is `A`
#'  - for `model = '1compartment'`, the parameters estimated are `Vdist`,
#'`kelim`, and possibly `kgutabs` and `Fgutabs` or `Fgutabs_Vdist` (see below).
#'  - for `model = '2compartment'`, the parameters estimated are `V1`, `kelim`,
#'`k12`, `k21`, and possibly `kgutabs` and `Fgutabs` or `Fgutabs_Vdist`.
#'
#'For 1-compartment and 2-compartment models, `kgutabs` will be estimated from
#'the data only if `fitdata` includes oral dosing data; otherwise it will be set
#'to NA.
#'
#'`Fgutabs` will be estimated from the data only if `fitdata` includes both oral
#'and IV data. If `fitdata` includes oral data but not IV data, then
#'`Fgutabs_Vdist` (1-compartment) or `Fgutabs_V1` (2-compartment) will be
#'estimated instead of `Fgutabs` and `Vdist` or `V1` (because only the ratio of
#'`Fgutabs` and `Vdist` is identifiable in that case).
#'
#'In addition to the model parameters, the standard deviation of the residual
#'errors will be estimated. If `fitdata` includes more than one unique study
#'value (`length(unique(fitdata$Study))>1`) and `pool_sigma == FALSE`, then a
#'separate error SD will be estimated for each unique study. These error SDs
#'will be named following the pattern `sigma_study_StudyID`, where `StudyID` is
#'one unique value of `fitdata$Study`, coerced to character if necessary. This
#'reflects an assumption that all concentration vs. time studies in `fitdata`
#'obey the same underlying PK model, but each study may have a different amount
#'of random measurement error (analogous to meta-regression with fixed effects).
#'
#'# Parameter estimation via numerical optimization
#'
#'Parameters are estimated using [optimx::optimx()].
#'
#'The objective function to be minimized is the negative log-likelihood, defined
#'in [log_likelihood()] with argument `negative = TRUE`. See documentation for
#'[log_likelihood()] for details.
#'
#'Parameter starting values are set in [get_starts()]. Parameter lower bounds
#'are set in [get_lower_bounds()]. Parameter upper bounds are set in
#'[get_upper_bounds()].
#'
#'By default, optimization is performed using the "bobyqa" method (see
#'[minqa::bobyqa()]). The method to use is specified in argument
#'`optimx_args$method`.
#'
#'Because parameter bounds are used, `optimx_args$method` should name a method
#'supported in [optimx::optimx()] that supports parameter bounds (i.e., one of
#'"bobyqa", "L-BFGS-B", "nlminb", "spg", "Rcgmin", "Rvmmin", "nmkb", "hjkb"). If
#'a different method is specified, then parameter bounds will be ignored, and
#'resulting parameter estimates may be non-physical.
#'
#'The objective function gradient is estimated numerically.
#'
#'## Parameter standard deviations via numerical Hessian
#'
#'Parameter standard deviations (uncertainty) are estimated from the Hessian
#'matrix (matrix of second derivatives of the objective function). The Hessian
#'matrix is estimated numerically. Parameter SDs are estimated as the square
#'root of the diagonal of the inverse of the Hessian (the matrix of second
#'derivatives of the objective function). The inverse Hessian approximates the
#'variance-covariance matrix of estimated parameters. The square root of its
#'diagonal approximates parameter standard deviations.
#'
#'If the numerically-estimated Hessian cannot be inverted, a pseudoinverse is
#'attempted.
#'
#'If the square root of the inverse Hessian contains NaN, then if
#'`optimx_args$method %in% 'L-BFGS-B'`, optimization is repeated with a smaller
#'convergence tolerance, until either there are no more NaNs or convergence
#'tolerance factor is equal to 1 (i.e. tolerance is at machine epsilon). This is
#'not done if `optimx_args$method %in% 'bobyqa'` because that method does not
#'use a convergence tolerance control parameter.
#'
#'# Expected variables in \code{fitdata}
#'
#' \describe{\item{\code{Time}}{Time for each data point}
#' \item{\code{Value}}{Concentration for each data point}
#' \item{\code{Dose}}{Dose for each data point} \item{\code{DTXSID}}{DSSTox
#' Substance ID identifying the substance being dosed. This function expects
#' only one unique value in this column; if there are multiple values, it will
#' stop with an error.} \item{\code{Species}}{String identifying the species
#' being dosed. This function expects only one unique value in this column; if
#' there are multiple values, it will stop with an error.}
#' \item{\code{Study}}{String identifying the study or study for each
#' data point. Each study is assumed to have its own residual error standard
#' deviation.} }
#'
#'@param fitdata A \code{data.frame} containing the set of concentration vs.
#'  time data to be fitted. See Details for expected variables. See also
#'  `preprocess_data()`.
#'@param model Which general model should be fit for each chemical. Presently,
#'  only "flat", "1compartment", and "2compartment" are implemented.
#'@param modelfun Either "analytic" or "full" -- whether to fit using the
#'  analytic solution to the model, or the full ODE model. Presently, "analytic"
#'  is recommended (because the analytic solution is exact and much faster). For
#'  `modelfun = 'full'`, only `model = '1compartment'` is supported.
#'@param pool_sigma Logical: Whether to pool all data (estimate only one error
#'  standard deviation) or not (estimate separate error standard deviations for
#'  each study). Default FALSE to estimate separate error SDs for each study.
#'  (If `fitdata` only includes one study, `pool_sigma` will have no effect,
#'  because only one error SD would be estimated in the first place.)
#'@param fit_conc_dose Logical: Whether to fit dose-normalized concentrations
#'  (TRUE) or non-dose-normalized concentrations (FALSE). Default TRUE.
#'@param rescale_time Logical: Whether to rescale time to make scale of time
#'  constants better-behaved. If TRUE, then if maximum time is more than 72
#'  hours, time will be converted to days before fitting; if maximum time is
#'  more than 1440 hours, time will be converted to months before fitting.
#'@param get_starts_args Any additional arguments to [get_starts()] (other than
#'  `model` and `fitdata`, which are always passed). Default NULL to accept the
#'  default arguments for [get_starts()].
#'@param get_lower_args Any additional arguments to [get_lower_bounds()] (other
#'  than `model` and `fitdata`, which are always passed). Default NULL to accept
#'  the default arguments for [get_lower_bounds()].
#'@param get_upper_args Any additional arguments to [get_upper_bounds()] (other
#'  than `model` and `fitdata`, which are always passed). Default NULL to accept
#'  the default arguments for [get_upper_bounds()].
#'@param optimx_args A named list of additional arguments to [optimx::optimx()],
#'  other than `par`, `fn`, `lower`, and `upper`. Default is:
#'
#'    ```
#'     list(
#'           "method" = "bobyqa",
#'           "itnmax" = 1e6,
#'          "control" = list("kkt" = FALSE))
#'    ```
#'  Briefly:  `"method"` allows you to select an optimization algorithm.
#'  `"itnmax"` controls the maximum number of iterations allowed in an attempt
#'  to optimize. `"control"` is itself a named list of control parameters
#'  relevant to the selected method. See documentation for [optimx::optimx()]
#'  for more details and more options. Note lower and upper bounds (box
#'  constraints) will be supplied; if you want them to be respected, please
#'  choose a method that allows box constraints (e.g. "bobyqa" or "L-BFGS-B").
#'
#'@return A `data.frame` of model-fitting results.
#'@author Caroline Ring, Chris Cook, John Wambaugh

analyze_subset <- function(fitdata,
                           model,
                           modelfun,
                           pool_sigma = FALSE,
                           fit_conc_dose = TRUE,
                           rescale_time = TRUE,
                           get_starts_args = list(start_from_httk = "all",
                                                  start_from_data = "all"),
                           get_lower_args = list(Vdist_from_species = FALSE),
                           get_upper_args = list(Fgutabs_Vdist_from_species = FALSE,
                                                 sigma_from_data = TRUE),
                           optimx_args = list(
                             "method" = "bobyqa",
                             "itnmax" = 1e6,
                             "control" = list("kkt" = FALSE)
                           ),
                           suppress.messages = FALSE){

  #convert back to data.frame
  fitdata <- as.data.frame(fitdata)

  this.dtxsid <- unique(fitdata$DTXSID)
  if(length(this.dtxsid) > 1) stop("analyze_subset(): More than one DTXSID in data")

  this.species <- unique(fitdata$Species)
  if(length(this.species) > 1) stop("analyze_subset(): More than one species in data")

  #check whether there is more than one study or not
  nstudy <- length(unique(fitdata$Study))
  studies_analyzed <- paste(sort(unique(fitdata$Study)),
                            collapse =", ")
  refs_analyzed <- paste(sort(unique(fitdata$Reference)),
                         collapse =", ")
  if (nstudy>1) {
    if(pool_sigma %in% FALSE){
      analysis_type <- "Joint Analysis"
    }else{
      analysis_type <- "Pooled Analysis"
    }
  }else{
    analysis_type <- "Single-Study Analysis"
  }

  n_subj <- range(fitdata$N_Subjects)
  if(length(unique(n_subj))==1) n_subj <- unique(n_subj)

  n_routes <- paste("iv: ",
                  sum(fitdata$Route %in% "iv"),
                  "; po: ",
                  sum(fitdata$Route %in% "po"))

  n_media <- paste("blood: ",
                   sum(fitdata$Media %in% "blood"),
                   "; plasma: ",
                   sum(fitdata$Media %in% "plasma"))

  if(!suppress.messages){
    message(paste0("Beginning ",
                   tolower(analysis_type),
    "for:\n",
                  this.dtxsid, " ", unique(fitdata$Compound), "\n",
                   "Species = ", this.species, "\n",
                  "Reference IDs = ",
                  refs_analyzed,
                  "\n",
                   "Study IDs = ",
                   studies_analyzed,
                   "\n",
                   "Number of observations = ",
                   nrow(fitdata),
                  " (iv: ",
                  sum(fitdata$Route %in% "iv"),
                  "; po: ",
                  sum(fitdata$Route %in% "po"),
                  ")\n",
                  "Number of subjects per observation = ",
                  paste(n_subj,
                        collapse = "-")
    )
    )
  }

  fitdata$Time.Hours <- fitdata$Time
  #rescale time if so requested
  if(rescale_time %in% TRUE){
    last_detect_time <- max(fitdata[is.finite(fitdata$Value),
                                    "Time.Hours"])

    if(last_detect_time > (24*365*2)){
      new_time_units <- "years"
    }else if(last_detect_time > (24*30*2)){
      new_time_units <- "months"
    }else if(last_detect_time > (24*7*2)){
      new_time_units <- "weeks"
    }else if(last_detect_time > 24*2){
      new_time_units <- "days"
    }else if(last_detect_time<0.5){
      new_time_units <- "minutes"
    }else{
      new_time_units = "hours"
    }

    #use lubridate to handle conversion of hours into new_time_units
    fitdata$Time <- convert_time(x = fitdata$Time.Hours,
                                 from = "hours",
                                 to = new_time_units,
                                 inverse = FALSE)

    fitdata$Time.Units <- new_time_units

    if(!(suppress.messages %in% TRUE) &
       !(new_time_units %in% "hours")){
      message(paste("Rescaling time from hours to",
                    new_time_units, "\n",
                    "Latest detection time = ",
                    signif(last_detect_time, 3),
                    "hours.",
                    "\nOld time range: ",
                    paste(signif(range(fitdata$Time.Hours), 3),
                          collapse = "-"),
                    "hours",
                    "\nNew time range:",
                    paste(signif(range(fitdata$Time), 3),
                          collapse = "-"),
                    new_time_units))
    }
  }

  #get parameter names and units, and determine whether to optimize each of
  #these parameters or not
  par_DF <- do.call(get_opt_params,
                    list("model" = model,
                         "fitdata" = fitdata,
                         "pool_sigma" = pool_sigma,
                         "suppress.messages" = suppress.messages))

  n_opt <- sum(par_DF$optimize_param %in% TRUE)

  #assign optimx::optimx() default max iterations/function evals if missing
  #this is default from optimx code itself
  if(is.null(optimx_args$itnmax)){
    optimx_args$itnmax <- max(5e4,
                      5000*round(
                        sqrt(
                          n_opt+1
                        )
                      )
        )


    if(!suppress.messages){
      message(paste("optimx_args does not include item 'itnmax'.",
                    "Setting optimx_args$itnmax =",
                    "max(5e4, 5000*round(sqrt(n+1)))",
                    "where n = number of params to be optimized =",
                    n_opt,
                    "(This is the optimx() default)"))
    }
  }

  #assign convergence tolerance factor if missing and relevant
  if(optimx_args$method %in% "L-BFGS-B" &
     is.null(optimx_args$control$factr)){
    optimx_args$control$factr <- 1e7
    if(!suppress.messages){
      message(paste("Optimization method is",
                    "L-BFGS-B",
                    "but optimx::optimx() control argument",
                    "'factr' was not provided;",
                    "setting optimx_args$control$factr = 1e7"))
    }
  }

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

   out_DF$time_units_fitted <- "hours"

   out_DF[, fitted_types] <- NA_real_

   out_DF$Studies.Analyzed <- studies_analyzed
   out_DF$References.Analyzed <- refs_analyzed
   out_DF$Data.Analyzed <- analysis_type

   #fill in the loglike and AIC with NA s since no fit was done
   out_DF$LogLikelihood <-  NA_real_
   out_DF$AIC <-  NA_real_


   out_DF$flag <- NA_character_
   #Record the unique routes in this dataset
   #Route info provides context for why some parameters were/were not estimated
   out_DF$N_Routes <- n_routes
   #Record the unique media in this dataset
   out_DF$N_Media <- n_media

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
                       "fit_conc_dose" = fit_conc_dose,
                       "suppress.messages" = suppress.messages),
                     get_upper_args)
 )

 #get starting values
 par_DF <- do.call(get_starts,
                   c(list("par_DF" = par_DF,
                       "model" = model,
                       "fitdata" = fitdata,
                       "pool_sigma" = pool_sigma,
                       "fit_conc_dose" = fit_conc_dose,
                       "suppress.messages" = suppress.messages),
                     get_starts_args))

 #Initialize out_DF
 #There will be one row for each parameter
out_DF <- par_DF[par_DF$optimize_param %in% TRUE, ]
#add a variable for time units
out_DF$time_units_fitted <- new_time_units

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
    "..."))

    message(paste("Initial values:    ",
                  paste(apply(data.frame(Names = names(opt_params),
                                         Values = unlist(opt_params),
                                         stringsAsFactors = F),
                              1, function(x) paste(x, collapse = ": ")),
                        collapse = ", "),
                  sep = ""))

    #if a convergence tolerance factor is relevant,
    #tell us what it is -- otherwise don't
    if(optimx_args$method %in% "L-BFGS-B"){
      message(paste("Convergence tolerance factor =",
                     optimx_args$control$factr))
    }
  }

  all_data_fit <- tryCatch({

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
          fit_conc_dose = fit_conc_dose,
          force_finite = (optimx_args$method %in% "L-BFGS-B"),
          negative = TRUE
        ) #end list()
      ) #end args = c()
    ) #end do.call

    #output of optimx::optimx():
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
  tmp <- try(any(is.na(as.vector(stats::coef(all_data_fit)))))
  if(inherits(tmp, "try-error")){
    fitfail <- TRUE
  }else{
    fitfail <- FALSE
  }

  if(fitfail %in% TRUE){
  out_DF[, fitted_types] <- NA_real_


    out_DF$Studies.Analyzed <- studies_analyzed
    out_DF$References.Analyzed <- refs_analyzed
    out_DF$Data.Analyzed <- analysis_type

  #fill in the loglike and AIC with NA s since no fit was done
  out_DF$LogLikelihood <-  NA_real_
  out_DF$AIC <-  NA_real_


  out_DF$flag <- NA_character_
  #Record the unique routes in this dataset
  #Route info provides context for why some parameters were/were not estimated
  out_DF$N_Routes <- n_routes
  #Record the unique media in this dataset
  out_DF$N_Media <- n_media

  #include a message about why no fit was done
  msg <- paste("Optimization failed. Error message from optimx():",
               all_data_fit)
  out_DF$message <- msg

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

#if fit was successful:
#post-processing of fit results
  #

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
    log_likelihood(x,
                        const_params = const_params,
                        DF = fitdata,
                        modelfun = modelfun,
                        model = model,
                   fit_conc_dose = fit_conc_dose,
                        force_finite = (optimx_args$method %in% "L-BFGS-B"),
                   negative = TRUE)
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
                       #pseudovariance matrix
                       #see http://gking.harvard.edu/files/help.pdf
                       suppressWarnings(tmp <- tryCatch(
                         diag(chol(MASS::ginv(numhess),
                                        pivot = TRUE)) ^ (1/2),
                         error = function(err){
                           if (!suppress.messages) {
                             message(paste0("Pseudovariance matrix failed,",
                                            " returning NAs"))
                           }
                           rep(NA_real_, nrow(numhess))
                         }
                       )
                       )
                       return(tmp)
                     })
  names(sds) <- names(means)



  #If any of the SDs are NaN,
  #repeat optimization with smaller convergence tolerance
  #-- only if method is "L-BFGS-B" -- other methods do not use this
  if(optimx_args$method %in% "L-BFGS-B"){
  while(any(is.nan(sds)) & optimx_args$control$factr > 1) {


    #then redo optimization

    #start from fitted values

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

    }


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
                                 fit_conc_dose = fit_conc_dose,
                                 force_finite = (optimx_args$method %in% "L-BFGS-B"),
                                 negative = TRUE
                               )
                             )
    )

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
      log_likelihood(x,
                          const_params = const_params,
                          DF = fitdata,
                          modelfun = modelfun,
                          model = model,
                     fit_conc_dose = fit_conc_dose,
                          force_finite = (optimx_args$method %in% "L-BFGS-B"),
                     negative = TRUE)
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
                      #pseudovariance matrix
                      #see http://gking.harvard.edu/files/help.pdf
                      suppressWarnings(tmp <- tryCatch(
                        diag(chol(MASS::ginv(numhess),
                                  pivot = TRUE)) ^ (1/2),
                        error = function(err){
                          if (!suppress.messages) {
                            message(paste0("Pseudovariance matrix failed,",
                                           " returning NAs"))
                          }
                          rep(NA_real_, nrow(numhess))
                        }
                      )
                      )
                      return(tmp)
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
  par_DF$time_units_fitted <- new_time_units
  out_DF <- merge(par_DF,
                  fit_DF,
                  by = "param_name",
                  all.y = TRUE,
                  all.x = FALSE)

    out_DF$Studies.Analyzed <- studies_analyzed
    out_DF$References.Analyzed <- refs_analyzed
    out_DF$Data.Analyzed <- analysis_type

  #Add log-likelihood and AIC values

  #convert negative log-likelihood (the value
  #returned by optimizer) to regular log-likelihood (multiply by -1)
  out_DF$LogLikelihood <- -1 * all_data_fit$value
  out_DF$AIC <- 2 * length(opt_params) - 2 * -1 * all_data_fit$value

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
  out_DF$N_Routes <- n_routes
  #Record the unique media in this dataset
  out_DF$N_Media <- n_media
  out_DF$message <- "Optimization successful."

#Keep information about optimization
  #number of function evals
  out_DF$fevals <- as.integer(all_data_fit$fevals)
  #convergence code
  out_DF$convcode <- as.integer(all_data_fit$convcode)
  #number of iterations
  out_DF$niter <- as.integer(all_data_fit$niter)
  #optimziation method used
  out_DF$method <- optimx_args$method
  #optimx control parameters
  #record control params
  out_DF[paste0("control_",
                names(optimx_args$control))] <- optimx_args$control


return(out_DF)


}

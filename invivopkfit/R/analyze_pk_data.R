#' Actually does the fitting
#'
#' Fits model parameters to concentration vs. time data for a given chemical
#'
#'
#' @param fitdata A data.table of concentration vs. time data for a given
#' chemical
#' @param this.dtxsid A DSSTox Substance Identifier for the chemical to be fitted
#' @param paramnames A list of names of the model parameters to be fitted
#' @param modelfun "analytic" to use the analytic model solution, "full" to use
#' the full ODE model
#' @param model The model to fit, either "1compartment" or "2compartment"
#' (other models not implemented for now)
#' @param this.reference placeholder
#' @return A single row of fitted parameter values (arithmetic means, geometric
#' means, modes, arithmetic standard deviations, and geometric standard
#' deviations)
#' @author Caroline Ring, John Wambaugh
analyze_pk_data <- function(fitdata,
                            this.dtxsid,
                            paramnames,
                            modelfun,
                            model,
                            this.reference = NULL,
                            suppress.messages = FALSE) {
  UPPERBOUNDARY <- 1e4

  #take a copy of input data table so it behaves as though passed by value
  fitdata <- data.table::copy(fitdata)

  # If there is no reference column, add on in:
  if (!is.null(this.reference)) {
    fitdata[, Reference := this.reference]
  }

  ### if the subset of data has LOQ values of 'NA', make the LOQ values = 0.45 * the minimum Value
  fitdata[, LOQ := as.numeric(LOQ)]
  fitdata[is.na(LOQ), LOQ := 0.45 * min(Value, na.rm = TRUE), by = .(Compound, Reference, Media)]

  if (!suppress.messages) cat(paste("Optimizing data for chemical ", this.dtxsid, "\n", sep =""))

  #Set a plausible upper bound for sigma:
  TRYSIGMA <- stats::median(fitdata$Value, na.rm = T)
  MAXSIGMA <- 100 * stats::median(fitdata$Value, na.rm = T)

  #Initialize output table with HTTK-predicted parameter values
  out.dt <- fitdata[1, paramnames, with = FALSE]

  #
  #
  # ADD OTHER COMPOUND IDENTIFIERS
  #
  #

  out.dt[, Compound := fitdata[1,"Compound"]]
  out.dt[, CAS := fitdata[1,"CAS"]]
  out.dt[, Media := fitdata[1,"Media"]]

  #if (this.dtxsid == "DTXSID0020442") browser()

  # Describe initial parameter values
  out.dt[, param.value.type := 'Predicted']

  these.params <- as.list(fitdata[1, paramnames, with = FALSE])

  #Add a different sigma value for each reference
  refs <- fitdata[,
                  unique(Reference)]
  # Initialize the standard deviations to 0.2
  these.params[sapply(refs,
                      function(x) paste('sigma2',
                                        x,
                                        sep = '.'))] <- rep(max(TRYSIGMA / 100, 0.1),
                                                            length(refs))

  #log-transform the model parameters
  these.params <- lapply(these.params, log)

  ### this doesn't seem to be necessary
  # #Construct the log-likelihood for the data on the current chemical
  # out <- log.likelihood(params=these.params,
  #                       DT=fitdata,
  #                       modelfun=modelfun,
  #                       model=model)

  #get list of params to optimize
  if (model == '1compartment') {
    opt.params <- these.params[names(these.params) %in%
                                 c("Vdist",
                                   "kelim")]

  } else if (model == '2compartment') {
    opt.params <- these.params[names(these.params) %in%
                                 c("V1",
                                   "Ralphatokelim",
                                   "Fbetaofalpha",
                                   "kelim")]

  } else if (model == 'flat') {
      opt.params <- these.params[names(these.params) %in% "A"]
      opt.params["A"] <- log(mean(fitdata$Value, na.rm = TRUE))
    }

  if (model != 'flat') {

  # Guess a good value for Vd and kelim:
  if ("iv" %in% fitdata$Route)
  {
    elim.data <- subset(fitdata, Route == "iv" & !is.na(Value))
    elim.data$Value <- elim.data$Value / elim.data$Dose
    elim.data <- subset(elim.data, Dose > 0)
    opt.params["kelim"] <- log(max(-stats::lm(log(Value) ~ Time,data=elim.data)$coefficients["Time"], 0.0001))
    if (is.na(opt.params["kelim"])) opt.params["kelim"] <- log(1)
    Vd.data <- subset(elim.data, Time == min(Time))
    if (model == '1compartment')
    {
      opt.params["Vdist"] <- log(1 / mean(Vd.data$Value))

      # foo <- function(elim.data, opt.params) {
      #
      #   Value <- cp_1comp(time = elim.data$Time,
      #                     params = list(kelim = opt.params$kelim,
      #                                   Vdist = opt.params$Vdist,
      #                                   Fgutabs = opt.params$Fgutabs,
      #                                   kgutabs = opt.params$kgutabs),
      #                     dose = elim.data$Dose,
      #                     iv.dose = elim.data$iv)
      #
      # }
      #
      # optimx(par = unlist(opt.params), fn = foo(opt.params = opt.params, elim.data = elim.data))

      # tryCatch({
      #
      # one_comp_form <- formula(log(Value) ~ log(cp_1comp(time = Time,
      #                                                       params = list(kelim = kelim,
      #                                                                     Vdist = Vdist,
      #                                                                     Fgutabs = 1,
      #                                                                     kgutabs = 2.18),
      #                                                       dose = Dose,
      #                                                       iv.dose = iv[[1]])))
      #
      #     ### one comp iv fit
      #     ### use calculations for kelim and Vdist as starting values
      # one_comp_fit <- nls2::nls2(formula = one_comp_form, data = elim.data,
      #                               start = list(kelim = exp(opt.params$kelim),
      #                                            Vdist = exp(opt.params$Vdist)))
      #
      #     ### edit opt.params values to nls2 model output
      # opt.params["kelim"] <- log(coef(one_comp_fit)["kelim"])
      # opt.params["Vdist"] <- log(coef(one_comp_fit)["Vdist"])
      #
      #   },
      #   error = function(cond) {
      #     message("nls2 pre-optimization failed for iv data. Using less refined parameter estimates as starting values in optimizer.")
      #     return(NA)
      #   })

    } else if (model == '2compartment') {

      opt.params["V1"] <- log(1 / mean(Vd.data$Value))
      if (length(unique(elim.data$Time)) > 3) {

        elim.data.test <- elim.data

        elim.data <- subset(elim.data,Time <= sort(unique(Time))[3]) ### added sort
        alpha <- -stats::lm(log(Value) ~ Time, data = elim.data)$coefficients["Time"]
        opt.params["Ralphatokelim"] <- log(max(alpha / exp(opt.params[["kelim"]]), 2))

        # tryCatch({
        #
        # two_comp_form <- formula(log(Value) ~ log(cp_2comp(time = Time,
        #                                                       params = list(kelim = kelim,
        #                                                                     V1 = V1,
        #                                                                     Fgutabs = 1,
        #                                                                     kgutabs = 2.18,
        #                                                                     Fbetaofalpha = Fbetaofalpha,
        #                                                                     Ralphatokelim = Ralphatokelim),
        #                                                       dose = Dose,
        #                                                       iv.dose = iv[[1]])))
        #
        # ### two comp iv fit
        # ### use calculations for kelim and Vdist as starting values
        # ### opt.params is log scale
        # two_comp_fit <- nls2::nls2(formula = two_comp_form, data = elim.data.test,
        #                               start = list(kelim = exp(opt.params$kelim),
        #                                            V1 = exp(opt.params$V1),
        #                                            Fbetaofalpha = exp(opt.params$Fbetaofalpha),
        #                                            Ralphatokelim = exp(opt.params$Ralphatokelim)),
        #                               algorithm = "port",
        #                               lower = c(0, 0, 0, 1))
        #
        #
        #
        #
        # ### edit opt.params values to nls2 model output
        # opt.params["kelim"] <- log(coef(two_comp_fit)["kelim"])
        # opt.params["V1"] <- log(coef(two_comp_fit)["V1"])
        # opt.params["Fbetaofalpha"] <- log(coef(two_comp_fit)["Fbetaofalpha"])
        # opt.params["Ralphatokelim"] <- log(coef(two_comp_fit)["Ralphatokelim"])
        #     }
        # ,
        #
        #     error = function(cond) {
        #       message("nls2 pre-optimization failed for iv data. Using less refined parameter estimates as starting values in optimizer.")
        #       return(NA)
        #     },
        #
        # finally = {
        #   message("yippee")
        # })

      }
      else alpha <- log(2)
      opt.params["Fbetaofalpha"] <- log(0.25)
    }
  } else if ("po" %in% fitdata$Route) {
    # browser()
    elim.data <- subset(fitdata, Route == "po" & !is.na(Value))
    elim.data$Value <- elim.data$Value/elim.data$Dose
    elim.data <- subset(elim.data,Dose > 0)
    max.time <- as.numeric(min(elim.data[Value == max(Value), "Time"]))
    elim.data <- subset(elim.data, Time >= max.time)
    opt.params["kelim"] <- log(max(-stats::lm(log(Value) ~ Time, data = elim.data)$coefficients["Time"], 0.0001))
    if (is.na(opt.params["kelim"])) opt.params["kelim"] <- log(1)
    Vd.data <- subset(elim.data, Time == min(Time))
    if (model =='1compartment') {

      opt.params["Vdist"] <- log(1 / mean(Vd.data$Value))

    } else {
      opt.params["V1"] <- log(1 / mean(Vd.data$Value))
      opt.params["Ralphatokelim"] <- log(2)
      opt.params["Fbetaofalpha"] <- log(0.25)
    }
  }

  if (is.na(opt.params[["kelim"]])) opt.params["kelim"] <- log(10 ^ -5) =
    opt.params["kelim"] <- max(opt.params[["kelim"]], log(10 ^ -5))


  #if data includes po, then also optimize oral fraction absorbed and gut
  #absorption rate
  if ("po" %in% fitdata$Route) {
    opt.params <- c(opt.params,
                    these.params["kgutabs"])

    # if (model == "1compartment") {
    #   # tryCatch(
    #   #   {
    #   one_comp_form <- formula(log(Value) ~ log(cp_1comp(time = Time,
    #                                                         params = list(kelim = kelim,
    #                                                                       Vdist = Vdist,
    #                                                                       Fgutabs = Fgutabs,
    #                                                                       kgutabs = kgutabs),
    #                                                         dose = Dose,
    #                                                         iv.dose = iv)))
    #
    #   ### one comp po fit
    #   ### use calculations for kelim and Vdist as starting values
    #   one_comp_fit_po <- nls2::nls2(formula = one_comp_form_po, data = elim.data,
    #                                 start = list(kelim = exp(opt.params$kelim),
    #                                              Vdist = exp(opt.params$Vdist),
    #                                              # Fgutabs = 1,
    #                                              kgutabs = exp(opt.params$kgutabs)
    #       ))
    #
    #   ### edit opt.params values to nls2 model output
    #   opt.params["kelim"] <- log(coef(one_comp_fit)["kelim"])
    #   opt.params["Vdist"] <- log(coef(one_comp_fit)["Vdist"])
    #   opt.params["kgutabs"] <- log(coef(one_comp_fit)["kgutabs"])
    #
    #   # },
    #   # error = function(cond) {
    #   #   message("nls2 pre-optimization failed for iv data. Using less refined parameter estimates as starting values in optimizer.")
    #   #   return(NA)
    #   # })
    # }

    # Need oral and iv data to get at Fgutabs:
    if ("iv" %in% fitdata$Route) {

      opt.params <- c(opt.params,
                      these.params["Fgutabs"])
      iv.data <- subset(fitdata, Route == "iv" & !is.na(Value))
      iv.data <- subset(iv.data, Time == min(Time))
      iv.data <- subset(iv.data, Dose == max(Dose))
      oral.data <- subset(fitdata,Route == "po" & !is.na(Value))
      oral.data <- subset(oral.data, Value == max(Value))[1, ]
      opt.params["Fgutabs"] <- log(oral.data$Value / mean(iv.data$Value) * mean(iv.data$Dose) / mean(oral.data$Dose))
      if (is.na(opt.params[["Fgutabs"]]) | opt.params[["Fgutabs"]] > 0) opt.params["Fgutabs"] <- log(0.5)

      # browser()
      # if (model == "1compartment") {
      #   # tryCatch(
      #   #   {
      #   test.fitdata <- subset(fitdata, !is.na(Value))
      #   one_comp_form_po <- formula(log(Value) ~ log(cp_1comp(time = Time,
      #                                                         params = list(kelim = kelim,
      #                                                                       Vdist = Vdist,
      #                                                                       Fgutabs = Fgutabs,
      #                                                                       kgutabs = kgutabs),
      #                                                         dose = Dose,
      #                                                         iv.dose = FALSE)))
      #
      #   ### one comp po fit
      #   ### use calculations for kelim and Vdist as starting values
      #   one_comp_fit_po <- nls2::nls2(formula = one_comp_form_po, data = test.fitdata,
      #                                 start = list(kelim = exp(opt.params$kelim),
      #                                              Vdist = exp(opt.params$Vdist),
      #                                              Fgutabs = exp(opt.params$Fgutabs),
      #                                              kgutabs = exp(opt.params$kgutabs)
      #                                 ))
      #
      #   ### edit opt.params values to nls2 model output
      #   opt.params["kelim"] <- log(coef(one_comp_fit_po)["kelim"])
      #   opt.params["Vdist"] <- log(coef(one_comp_fit_po)["Vdist"])
      #   opt.params["Fgutabs"] <- log(coef(one_comp_fit_po)["Fgutabs"])
      #   opt.params["kgutabs"] <- log(coef(one_comp_fit_po)["kgutabs"])
      #
      #   # },
      #   # error = function(cond) {
      #   #   message("nls2 pre-optimization failed for iv data. Using less refined parameter estimates as starting values in optimizer.")
      #   #   return(NA)
      #   # })
      # }
    }
  }
}


  #Add the per-study standard deviation to the set of params to optimize
  opt.params <- c(opt.params,
                  these.params[regexpr("sigma2",
                                       names(these.params)) != -1])

  #If the number of parameters is >= the number of data points,
  #then throw back everything NA with a message,
  #because there is no point wasting time trying to fit them.
  if (length(opt.params) >= nrow(fitdata)) {
    tmp.out <- data.table(param.value.type = c("Fitted arithmetic mean",
                                               "Fitted arithmetic std dev",
                                               "Fitted geometric mean",
                                               "Fitted geometric std dev",
                                               "Fitted mode"))
    out.dt <- rbind(out.dt, tmp.out, fill = TRUE)
    out.dt[, sigma := as.character(NA)]
    out.dt[, value := as.double(NA)]

    if (is.null(this.reference)) {

      out.dt[, Reference:=paste(sort(unique(fitdata$Reference)),
                                sep = ", ",
                                collapse =", ")]
      out.dt[, Data.Analyzed:=Reference]
      out.dt[regexpr(",", Data.Analyzed) != -1, Data.Analyzed := 'Joint Analysis']
    } else {
      # out.dt[,Reference:=this.reference]
      out.dt[, Data.Analyzed := this.reference]
    }

    out.dt[, LogLikelihood := as.numeric(NA)]
    out.dt[, AIC := as.double(NA)]

    cat(paste("For chemical ", this.dtxsid, " there were ", length(opt.params),
              " parameters and only ",
              nrow(fitdata), " data points. Optimization aborted.\n", sep = ""))
    # browser()
    return(out.dt)
  }

  #get list of params to hold constant and not optimize
  const.params <- these.params[!(names(these.params) %in%
                                   names(opt.params))]

  #
  #
  #
  # THIS IS WHERE OPTIMIZER UPPER BOUNDS ARE SET
  #
  #
  #

  # if (unique(fitdata$Compound) == "1,4-dioxane") browser()
  #change from a named list of params to optimize to a named vector of params to
  #optimize
  upper <- unlist(opt.params)

  #specify upper bounds of params to optimize (on a log scale!!)
  upper[] <- log(UPPERBOUNDARY)
  upper[regexpr("sigma",names(upper)) != -1] <- log(TRYSIGMA) ### was MAXSIGMA
  if (model == "1compartment") {
    upper["Vdist"] <- log(max(fitdata$Dose) / min(fitdata$LOQ))
  } else if (model == '2compartment') {
    upper["V1"] <- log(max(fitdata$Dose) / min(fitdata$LOQ))
    upper["Ralphatokelim"] <- log(1000)
    upper["Fbetaofalpha"] <- log(0.75) #on a log scale!
  }
  if ("Fgutabs" %in% unlist(opt.params)) {
    upper["Fgutabs"] <- log(1) #on a log scale!
  }
  if ("kgutabs" %in% unlist(opt.params)) {
    upper["kgutabs"] <- log(1000)
  }
# upper <- upper["Vdist"]
  # Force initial values to be within bounds:
  opt.params[opt.params>upper] <- upper[opt.params>upper] - 0.1

  #
  #
  #
  # THIS IS WHERE OPTIMIZER LOWER BOUNDS ARE SET
  #
  #
  #

  # Default lower bound of 10^-8 except for parameters where this wouldn't make sense:
  lower <- rep(log(10 ^ -8), length(upper))
  names(lower) <- names(upper)
  if (model == '1compartment') {
    lower["Vdist"] <- log(0.01)
  } else if (model == '2compartment') {
    lower["V1"] <- log(0.01)
    lower["Ralphatokelim"] <- 0
  }
  if ("Fgutabs" %in% names(lower)) {
    lower["Fgutabs"] <- log(0.05)
  }

  # Curve fitting doesn't really work if we let the standard deviation of the measurements get too small:
  lower[regexpr("sigma", names(lower)) != -1] <- log(0.00001)
  # Force initial values to be within bounds:
  opt.params[opt.params < lower] <- lower[opt.params < lower] + 0.1

  orig.params <- opt.params

  #factr: controls the convergence of the "L-BFGS-B" method. Convergence occurs when
  #the reduction in the objective is within this factor of the machine tolerance
  #(from ?optim)
  factr <- 1e7

  #Set up objective function to minimize:
  #minimize negative log-likelihood = maximize log-likelihood
  objfun <- function(x) {
    foo <- -log_likelihood(params = c(x,const.params),
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

  if (!suppress.messages) cat(paste("Initial values:    ", paste(apply(data.frame(Names = names(lapply(opt.params, exp)),
                                                          Values = unlist(lapply(opt.params, exp)),
                                                          stringsAsFactors = F),
                                               1, function(x) paste(x, collapse=": ")),
                                         collapse = ", "), "\n", sep = ""))
# browser()
  tryCatch(all.data.fit <- optimx::optimx(par = unlist(opt.params),
                                          fn = objfun,
                                          # lower=lower,
                                          upper=upper,
                                          method = "L-BFGS-B", #box constraints
                                          hessian = FALSE,
                                          control = list(factr = factr)),
           error = function(err){
             browser() #kick to debugger to find out what went wrong
           })

  #Get MLE params
  #ln.means <- all.data.fit$par
  ln.means <- as.vector(stats::coef(all.data.fit))
  names(ln.means) <- names(opt.params)

  # if (any(ln.means>upper)) browser()
  if (!suppress.messages) cat(paste("Optimized values:  ", paste(apply(data.frame(Names = names(ln.means),
                                                          Values = sapply(ln.means,exp),
                                                          stringsAsFactors = F),
                                               1, function(x) paste(x, collapse = ": ")),
                                         collapse = ", "),"\n", sep = ""))
  #Get SDs from Hessian
  #Calculate Hessian using function from numDeriv
  # browser()
  numhess <- numDeriv::hessian(func = objfun,
                               x = ln.means,
                               method = 'Richardson')
  ln.sds <- tryCatch(diag(solve(numhess)) ^ (1/2),
                     error = function(err){
                       #if hessian can't be inverted
                       if (!suppress.messages) cat("Hessian can't be inverted, using pseudovariance matrix to estimate parameter uncertainty.\n")
                       return(diag(chol(MASS::ginv(numhess),
                                        pivot = TRUE)) ^ (1/2)) #pseduovariance matrix
                       #see http://gking.harvard.edu/files/help.pdf
                     })
  names(ln.sds) <- names(ln.means)

  while(any(is.nan(ln.sds)) & factr > 1) {#If any of the SDs are NaN

    if (!suppress.messages) cat("One or more parameters has NaN standard deviation, repeating optimization with smaller convergence tolerance.\n")
    #then redo optimization
    #  if ("Fgutabs" %in% names(opt.params)) browser()
    opt.params <- ln.means + stats::runif(length(ln.means), -0.1, 0.1)
    opt.params[regexpr("sigma", names(opt.params)) != -1] <- opt.params[regexpr("sigma", names(opt.params))!=-1]+1
    opt.params[opt.params <= lower] <- lower[opt.params <= lower] + 0.1
    opt.params[opt.params >= upper] <- upper[opt.params >= upper] - 0.1

    if (!suppress.messages) cat(paste("Initial values:    ", paste(apply(data.frame(Names = names(opt.params),
                                                           Values = unlist(lapply(opt.params,exp)),
                                                           stringsAsFactors = F),
                                                1, function(x) paste(x, collapse = ": ")),
                                          collapse = ", "), "\n", sep = ""))
    factr <- factr / 10 #reducing factr by 10 (requiring closer convergence)

    #use general-purpose optimizer to optimize 1-compartment params to fit data
    #optimize by maximizing log-likelihood
    all.data.fit <- optimx::optimx(par = unlist(opt.params),
                                 fn = objfun,
                                 # lower=lower,
                                 upper=upper,
                                 method = "L-BFGS-B",
                                 hessian = FALSE,
                                 control = list(factr = factr))

    ln.means <- as.vector(stats::coef(all.data.fit))
    names(ln.means) <- names(opt.params)
    if (!suppress.messages) cat(paste("Optimized values:  ",paste(apply(data.frame(Names = names(ln.means),
                                                           Values = sapply(ln.means, exp),
                                                           stringsAsFactors=F),
                                                1, function(x) paste(x, collapse=": ")),
                                          collapse = ", "), "\n", sep = ""))

    #Get SDs from Hessian
    #Calculate Hessian using function from numDeriv
    numhess <- numDeriv::hessian(func = objfun,
                                 x = ln.means,
                                 method = 'Richardson')
    ln.sds <- tryCatch(diag(solve(numhess)) ^ (1/2),
                       error = function(err){
                         #if hessian can't be inverted
                         if (!suppress.messages) cat("Hessian can't be inverted, using pseudovariance matrix to estimate parameter uncertainty.\n")
                         return(diag(chol(MASS::ginv(numhess),
                                          pivot = TRUE)) ^ (1/2)) #pseduovariance matrix
                         #see http://gking.harvard.edu/files/help.pdf
                       })
    names(ln.sds) <- names(ln.means)
  }

  ### if(factr > 1) { this is to avoid unnecessary computation, but have not figured out where to put loop yet.
  # browser()
  # does removing this keep sigma in output?
  #Exclude sigma2 from results
  sigmas <- ln.means[regexpr("sigma", names(ln.means)) != -1]
  # ln.means <- ln.means[regexpr("sigma",names(ln.means))==-1]
  # ln.sds <- ln.sds[regexpr("sigma",names(ln.sds))==-1]

  # Actual means (tends to explode with large SD):
  #"Actual means" --> "arithmetic means"
  #as opposed to geometric means
  #the actual, non-log-transformed parameters are log-normally distributed
  #(the log-transformed parameters are normally distributed)
  #exp(ln.means) = geometric mean
  #exp(ln.sds) = geometric sd
  #have to be converted into arithmetic mean and sd.
  arith.means <- as.data.table(as.list(exp(ln.means +
                                             ln.sds ^ 2 / 2)))
  names(arith.means) <- names(ln.means)
  arith.means[, param.value.type := 'Fitted arithmetic mean']
  out.dt <- rbind(out.dt, arith.means, fill = TRUE)
  #If Fgutabs and kgutabs weren't optimized (because iv data),
  #then add them as NA to the output list.


  arith.se <- as.data.table(as.list(((exp(ln.sds ^ 2) - 1)*
                                       exp(2 * ln.means+
                                             ln.sds ^ 2)) ^ (1/2)))
  names(arith.se) <- names(ln.means)
  arith.se[, param.value.type := 'Fitted arithmetic std dev']
  out.dt <- rbind(out.dt, arith.se, fill = TRUE)
  # Geometric means and se:
  geo.means <- as.data.table(as.list(exp(ln.means)))
  names(geo.means) <- names(ln.means)
  geo.means[, param.value.type := 'Fitted geometric mean']
  out.dt <- rbind(out.dt, geo.means, fill = TRUE)

  geo.se <- as.data.table(as.list(ln.sds))
  names(geo.se) <- names(ln.means)
  geo.se[, param.value.type := 'Fitted geometric std dev']
  out.dt <- rbind(out.dt, geo.se, fill = TRUE)

  #Mode (the "peak" of the log-normal distribution)
  modes <- as.data.table(as.list(exp(ln.means - ln.sds ^ 2)))
  names(modes) <- names(ln.means)
  modes[, param.value.type := 'Fitted mode']
  out.dt <- rbind(out.dt, modes, fill = TRUE)

  out.dt <- as.data.table(tidyr::pivot_longer(out.dt,
                                              cols = grep("sigma", names(out.dt)),
                                              names_to = "sigma_id",
                                              values_to = "sigma_value"))

  if (is.null(this.reference)) {
    out.dt[, Reference := paste(sort(unique(fitdata$Reference)),
                              sep = ", ",
                              collapse = ", ")]
    out.dt[, Data.Analyzed := Reference]
    out.dt[regexpr(",", Data.Analyzed) != -1, Data.Analyzed := 'Joint Analysis']
  } else out.dt[, Data.Analyzed := this.reference]

  # If any of the parameters were not optimized or if the the model does not fit well:
  #  if (any(ln.means==as.vector(opt.params[names(ln.means)])) | any(sigmas>1))
  if (any(ln.means == as.vector(orig.params[names(ln.means)])) |
      any(sigmas > MAXSIGMA) |
      factr == 1 |
      all(ln.sds) == 0) { ### added factr and ln.sds

    out.dt[, LogLikelihood := 0]  ### replace with NA to make error more apparent
    out.dt[, AIC := Inf]
    if (!suppress.messages) cat(
      paste("Some parameters were not optimized, sigma >",
            MAXSIGMA,
            " some parameters had NaN standard deviation under smallest convergence tolerance, or the standard deviation of each parameter could not be calculated. Returning AIC=INF.\n"))
  } else if (model == '1compartment') {
    # Check for bad one-compartment model fits:
    if (ln.means["Vdist"] == upper["Vdist"]) {
      out.dt[ ,LogLikelihood := 0] ### replace with NA to make error more apparent
      out.dt[ ,AIC := Inf]
      if (!suppress.messages) cat("Vdist equal to upper bound on optimizer. Returning AIC=INF.\n")
    } else {
      out.dt[ ,LogLikelihood := -all.data.fit$value]
      out.dt[ ,AIC := 2 * length(opt.params) + 2 * all.data.fit$value]
    }
  } else if (model == '2compartment') {
    # Check for bad two-compartment model fits:
    if (ln.means["Fbetaofalpha"] == upper["Fbetaofalpha"] |
        ln.means["Fbetaofalpha"] == lower["Fbetaofalpha"] |
        ln.means["V1"] == upper["V1"])
    {
      out.dt[ ,LogLikelihood := 0] ### replace with NA to make error more apparent
      out.dt[ ,AIC := Inf]
      if (!suppress.messages) cat("Problem with Fbetaofalpha or V1 equaling optimizer bound. Returning AIC=INF.\n")
    } else {
      out.dt[, LogLikelihood := -all.data.fit$value]
      out.dt[, AIC := 2 * length(opt.params) + 2 * all.data.fit$value]
    }
  } else if (model == 'flat') {
    out.dt[, LogLikelihood := -all.data.fit$value]
    out.dt[, AIC := 2 * length(opt.params) + 2 * all.data.fit$value]
    out.dt[, Vdist := V1*(k12+k21)/k2]
  }

  # } else {warning("didn't converge")}

  if (!is.null(this.reference)) {
    out.dt[, Reference := NULL]
  }

  # Fill in missing Compound, CAS, media:
  out.dt[is.na(Compound), Compound := fitdata[1,"Compound"]]
  out.dt[is.na(CAS), CAS := fitdata[1,"CAS"]]
  out.dt[is.na(Media), Media := fitdata[1,"Media"]]

  #browser()
  return(out.dt)

}

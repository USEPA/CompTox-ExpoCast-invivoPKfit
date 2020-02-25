#'Actually does the fitting
#'
#'Fits model parameters to concentration vs. time data for a given chemical
#'
#'@param fitdata A data.table of concentration vs. time data for a given
#'  chemical
#'@param this.cas A CAS number for the chemical to be fitted
#'@param paramnames A list of names of the model parameters to be fitted
#'@param modelfun "analytic" to use the analytic model solution, "full" to use
#'  the full ODE model
#'@param model The model to fit, either "1compartment" or "2compartment" (other
#'  models not implemented for now)
#'
#' @author Caroline Ring, John Wambaugh
#'
#'@return A single row of fitted parameter values (arithmetic means, geometric
#'  means, modes, arithmetic standard deviations, and geometric standard
#'  deviations)

analyze.pk.data <- function(fitdata,
                            this.cas,
                            paramnames,
                            modelfun,
                            model,
                            this.reference=NULL) 
{
  UPPERBOUNDARY <- 1e4
  
  #take a copy of input data table so it behaves as though passed by value
  fitdata <- data.table::copy(fitdata)

  cat(paste("Optimizing data for CAS-RN ",this.cas,"\n",sep=""))

  #Set a plausible upper bound for sigma:
  MAXSIGMA <- 3*median(fitdata$Value,na.rm=T)

  #Initialize output table with HTTK-predicted parameter values
  out.dt <- fitdata[1,paramnames, with=FALSE]
  out.dt[, param.value.type:='Predicted']

  these.params <- as.list(fitdata[1, paramnames, with=FALSE])

  # If there is no reference column, add on in:
  if (!is.null(this.reference))
  {
    fitdata[,Reference:=this.reference]
  }

  #Add a different sigma value for each reference
  refs <- fitdata[,
                  unique(Reference)]
  # Initialize the standard devications to 0.2
  these.params[sapply(refs,
                      function(x) paste('sigma2',
                                        x,
                                        sep='.'))] <- rep(max(MAXSIGMA/100,
                                                            0.1),
                                                        length(refs))
  
  #log-transform the model parameters
  these.params <- lapply(these.params,log)

  #Construct the log-likelihood for the data on the current chemical
  out <- log.likelihood(params=these.params,
                        DT=fitdata,
                        modelfun=modelfun,
                        model=model)

  #get list of params to optimize
  if (model=='1compartment'){
    opt.params <- these.params[names(these.params) %in%
                                 c("Vdist",
                                   "kelim")]

  }else if (model=='2compartment'){
    opt.params <- these.params[names(these.params) %in%
                                 c("V1",
                                   "Ralphatokelim",
                                   "Fbetaofalpha",
                                   "kelim")]
  }

  # Guess a good value for Vd and kelim:
  if ("iv" %in% fitdata$Route)
  {
    elim.data <- subset(fitdata,Route=="iv"&!is.na(Value))
    elim.data$Value <- elim.data$Value/elim.data$Dose
    elim.data <- subset(elim.data,Dose>0)
    opt.params["kelim"] <- log(max(-lm(log(Value) ~ Time,data=elim.data)$coefficients["Time"],0.0001))
    if (is.na(opt.params["kelim"])) opt.params["kelim"] <- log(1)
    Vd.data <- subset(elim.data,Time==min(Time))
    if (model=='1compartment') 
    {
      opt.params["Vdist"] <- log(1/mean(Vd.data$Value))
    } else {
      opt.params["V1"] <- log(1/mean(Vd.data$Value))
      if (length(unique(elim.data$Time))>3)
      { 
        elim.data <- subset(elim.data,Time<=unique(Time)[3])
        alpha <- -lm(log(Value) ~ Time,data=elim.data)$coefficients["Time"]
        opt.params["Ralphatokelim"] <- log(max(alpha/exp(opt.params[["kelim"]]),2))
      }
      else alpha <- log(2)
      opt.params["Fbetaofalpha"] <- log(0.25)
    }
  } else if ("po" %in% fitdata$Route) {
    elim.data <- subset(fitdata,Route=="po"&!is.na(Value))
    elim.data$Value <- elim.data$Value/elim.data$Dose
    elim.data <- subset(elim.data,Dose>0)
    max.time <- as.numeric(min(elim.data[Value == max(Value), "Time"]))
    elim.data <- subset(elim.data,Time>=max.time)
    opt.params["kelim"] <- log(max(-lm(log(Value) ~ Time,data=elim.data)$coefficients["Time"],0.0001))
    if (is.na(opt.params["kelim"])) opt.params["kelim"] <- log(1)
    Vd.data <- subset(elim.data,Time==min(Time))
    if (model=='1compartment') 
    {
      opt.params["Vdist"] <- log(1/mean(Vd.data$Value))
    } else {
      opt.params["V1"] <- log(1/mean(Vd.data$Value))
      opt.params["Ralphatokelim"] <- log(2)
      opt.params["Fbetaofalpha"] <- log(0.25)
    }
  }
  
  if (is.na(opt.params[["kelim"]])) opt.params["kelim"]<-log(10^-5)
  opt.params["kelim"] <- max(opt.params[["kelim"]],log(10^-5))

  #if data includes po, then also optimize oral fraction absorbed and gut
  #absorption rate
  if ("po" %in% fitdata$Route)
  {
    opt.params <- c(opt.params,
                    these.params["kgutabs"])
    # Need oral and iv data to get at Fgutabs:
    if ("iv" %in% fitdata$Route)
    { 
      opt.params <- c(opt.params,
                      these.params["Fgutabs"])
      iv.data <- subset(fitdata,Route=="iv"&!is.na(Value))
      iv.data <- subset(iv.data,Time==min(Time))
      iv.data <- subset(iv.data,Dose==max(Dose))               
      oral.data <- subset(fitdata,Route=="po"&!is.na(Value))
      oral.data <- subset(oral.data,Value==max(Value))[1,]
      opt.params["Fgutabs"] <- log(oral.data$Value/mean(iv.data$Value)*mean(iv.data$Dose)/mean(oral.data$Dose))
      if (is.na(opt.params[["Fgutabs"]]) | opt.params[["Fgutabs"]]>0) opt.params["Fgutabs"]<-log(0.5)
    } 
  }

  
  #Add the per-study standard deviation to the set of params to optimize
  opt.params <- c(opt.params,
                  these.params[regexpr("sigma2",
                                       names(these.params))!=-1])

  #If the number of parameters is >= the number of data points,
  #then throw back everything NA with a message,
  #because there is no point wasting time trying to fit them.
  if (length(opt.params)>=nrow(fitdata)){
    tmp.out <- data.table(param.value.type=c("Fitted arithmetic mean",
                                             "Fitted arithmetic std dev",
                                             "Fitted geometric mean",
                                             "Fitted geometric std dev",
                                             "Fitted mode"))
    out.dt <- rbind(out.dt, tmp.out, fill=TRUE)

    if (is.null(this.reference))
    {
      out.dt[, Reference:=paste(sort(unique(fitdata$Reference)),
                                sep=", ",
                                collapse=", ")]
      out.dt[, Data.Analyzed:=Reference]
      out.dt[regexpr(",",Data.Analyzed)!=-1, Data.Analyzed:='Joint Analysis']
    } else {
      out.dt[, Data.Analyzed:=this.reference]
      out.dt[,Reference:=this.reference]
    }    
    out.dt[,LogLikelihood:=as.numeric(NA)]
    out.dt[,AIC:=as.numeric(NA)]
    cat(paste("For CAS ", this.cas, " there were ", length(opt.params),
                  " parameters and only ",
                  nrow(fitdata), " data points. Optimization aborted.\n",sep=""))

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
  
  #change from a named list of params to optimize to a named vector of params to
  #optimize
  upper <- unlist(opt.params)

  #specify upper bounds of params to optimize (on a log scale!!)
  upper[] <- log(UPPERBOUNDARY)
  upper[regexpr("sigma",names(upper))!=-1]<-log(MAXSIGMA) 
  if (model=='2compartment'){
    upper["Ralphatokelim"] <- log(1000)
    upper["Fbetaofalpha"] <- log(0.75) #on a log scale!
   }
  if ("Fgutabs" %in% unlist(opt.params)){
    upper["Fgutabs"] <- log(1) #on a log scale!
  }
  if ("kgutabs" %in% unlist(opt.params)){
    upper["kgutabs"] <- log(1000)
  }
  # Force initial values to be within bounds:
  opt.params[opt.params>upper] <- upper[opt.params>upper]-0.1
 
  #
  #
  #
  # THIS IS WHERE OPTIMIZER LOWER BOUNDS ARE SET
  #
  #
  #
  
  # Default lower bound of 10^-8 except for parameters where this wouldn't make sense:
  lower <- rep(log(10^-8), length(upper))
  names(lower) <- names(upper)
  if (model=='1compartment'){
    lower["Vdist"] <- log(0.01)
  } else if (model=='2compartment'){
    lower["V1"] <- log(0.01)
    lower["Ralphatokelim"] <- 0
  }
  if ("Fgutabs" %in% names(lower))
  {
    lower["Fgutabs"] <- log(0.05)
  }
  # Curve fitting doesn't really work if we let the standard deviation of the measurements get too small:
  lower[regexpr("sigma",names(lower))!=-1]<-log(0.00001)
  # Force initial values to be within bounds:
  opt.params[opt.params<lower] <- lower[opt.params<lower]+0.1
  
  orig.params <- opt.params

  #factr: controls the convergence of the "L-BFGS-B" method. Convergence occurs when
  #the reduction in the objective is within this factor of the machine tolerance
  #(from ?optim)
  factr <- 1e7

  #Set up objective function to minimize:
  #minimize negative log-likelihood = maximize log-likelihood
  objfun <- function(x) {
    foo <- -log.likelihood(params=c(x,const.params),
                                        DT=fitdata,
                                        modelfun=modelfun,
                                        model=model)
    #method L-BFGS-B requires finite values be returned,
    #so if log-likelihood is NA or -Inf,
    #just return a large negative log.likelihood
    #(= a large positive -log.likelihood)
    if (!is.finite(foo)) foo <- 99999
    return(foo)
  }

  cat(paste("Initial values:    ",paste(apply(data.frame(Names=names(lapply(opt.params,exp)),Values=unlist(lapply(opt.params,exp)),stringsAsFactors=F),1,function(x) paste(x,collapse=": ")),collapse=", "),"\n",sep=""))
  
  tryCatch(all.data.fit <- optimx::optimx(par=unlist(opt.params),
                                          fn=objfun,
                                          lower=lower,
                                          upper=upper,
                                          method="L-BFGS-B", #box constraints
                                          hessian = FALSE,
                                          control=list(factr=factr)),
           error = function(err){
             browser() #kick to debugger to find out what went wrong
           })

  #Get MLE params
  #ln.means <- all.data.fit$par
  ln.means <- as.vector(coef(all.data.fit))
  names(ln.means) <- names(opt.params)

  if (any(ln.means>upper)) browser()
  cat(paste("Optimized values:  ",paste(apply(data.frame(Names=names(ln.means),Values=sapply(ln.means,exp),stringsAsFactors=F),1,function(x) paste(x,collapse=": ")),collapse=", "),"\n",sep=""))
  #Get SDs from Hessian
  #Calculate Hessian using function from numDeriv
  numhess <- numDeriv::hessian(func=objfun,
                               x=ln.means,
                               method='Richardson')
  ln.sds <- tryCatch(diag(solve(numhess))^(1/2),
                     error = function(err){
                       #if hessian can't be inverted
                       cat("Hessian can't be inverted, using pseudovariance matrix to estimate parameter uncertainty.\n")
                       return(diag(chol(MASS::ginv(numhess),
                                        pivot=TRUE))^(1/2)) #pseduovariance matrix
                       #see http://gking.harvard.edu/files/help.pdf
                     })
  names(ln.sds) <- names(ln.means)

  while(any(is.nan(ln.sds)) & factr>1) #If any of the SDs are NaN
  {
    cat("One or more paramters has NaN standard deviation, repeating optimization with smaller convergence tolerance.\n")
    #then redo optimization
  #  if ("Fgutabs" %in% names(opt.params)) browser()
    opt.params <- ln.means+runif(length(ln.means),-0.1,0.1)
    opt.params[regexpr("sigma",names(opt.params))!=-1] <- opt.params[regexpr("sigma",names(opt.params))!=-1]+1
    opt.params[opt.params<=lower] <- lower[opt.params<=lower]+0.1
    opt.params[opt.params>=upper] <- upper[opt.params>=upper]-0.1
    
    cat(paste("Initial values:    ",paste(apply(data.frame(Names=names(opt.params),Values=unlist(lapply(opt.params,exp)),stringsAsFactors=F),1,function(x) paste(x,collapse=": ")),collapse=", "),"\n",sep=""))
    factr <- factr/10 #reducing factr by 10 (requiring closer convergence)

    #use general-purpose optimizer to optimize 1-compartment params to fit data
    #optimize by maximizing log-likelihood
    all.data.fit<-optimx::optimx(par=unlist(opt.params),
                                 fn=objfun,
                                 lower=lower,
                                 upper=upper,
                                 method="L-BFGS-B",
                                 hessian = FALSE,
                                 control=list(factr=factr)) 

    ln.means <- as.vector(coef(all.data.fit))
    names(ln.means) <- names(opt.params)
    cat(paste("Optimized values:  ",paste(apply(data.frame(Names=names(ln.means),Values=sapply(ln.means,exp),stringsAsFactors=F),1,function(x) paste(x,collapse=": ")),collapse=", "),"\n",sep=""))

    #Get SDs from Hessian
    #Calculate Hessian using function from numDeriv
    numhess <- numDeriv::hessian(func=objfun,
                                 x=ln.means,
                                 method='Richardson')
    ln.sds <- tryCatch(diag(solve(numhess))^(1/2),
                       error = function(err){
                         #if hessian can't be inverted
                         cat("Hessian can't be inverted, using pseudovariance matrix to estimate parameter uncertainty.\n")
                         return(diag(chol(MASS::ginv(numhess),
                                          pivot=TRUE))^(1/2)) #pseduovariance matrix
                         #see http://gking.harvard.edu/files/help.pdf
                       })
    names(ln.sds) <- names(ln.means)
  }

  #Exclude sigma2 from results
  sigmas <- ln.means[regexpr("sigma",names(ln.means))!=-1] 
  ln.means <- ln.means[regexpr("sigma",names(ln.means))==-1]
  ln.sds <- ln.sds[regexpr("sigma",names(ln.sds))==-1]

  # Actual means (tends to explode with large SD):
  #"Actual means" --> "arithmetic means"
  #as opposed to geometric means
  #the actual, non-log-transformed parameters are log-normally distributed
  #(the log-transformed parameters are normally distributed)
  #exp(ln.means) = geometric mean
  #exp(ln.sds) = geometric sd
  #have to be converted into arithmetic mean and sd.
  arith.means <- as.data.table(as.list(exp(ln.means +
                                             ln.sds^2/2)))
  names(arith.means) <- names(ln.means)
  arith.means[, param.value.type:='Fitted arithmetic mean']
  out.dt <- rbind(out.dt, arith.means, fill=TRUE)
  #If Fgutabs and kgutabs weren't optimized (because iv data),
  #then add them as NA to the output list.


  arith.se <- as.data.table(as.list(((exp(ln.sds^2) - 1)*
                                       exp(2*ln.means+
                                             ln.sds^2))^(1/2)))
  names(arith.se) <- names(ln.means)
  arith.se[, param.value.type:='Fitted arithmetic std dev']
  out.dt <- rbind(out.dt, arith.se, fill=TRUE)
  # Geometric means and se:
  geo.means <- as.data.table(as.list(exp(ln.means)))
  names(geo.means) <- names(ln.means)
  geo.means[, param.value.type:='Fitted geometric mean']
  out.dt <- rbind(out.dt, geo.means, fill=TRUE)

  geo.se <- as.data.table(as.list(ln.sds))
  names(geo.se) <- names(ln.means)
  geo.se[, param.value.type:='Fitted geometric std dev']
  out.dt <- rbind(out.dt, geo.se, fill=TRUE)

  #Mode (the "peak" of the log-normal distribution)
  modes <- as.data.table(as.list(exp(ln.means - ln.sds^2)))
  names(modes) <- names(ln.means)
  modes[, param.value.type:='Fitted mode']
  out.dt <- rbind(out.dt, modes, fill=TRUE)


  if (is.null(this.reference))
  {
    out.dt[, Reference:=paste(sort(unique(fitdata$Reference)),
                              sep=", ",
                              collapse=", ")]
    out.dt[, Data.Analyzed:=Reference]
    out.dt[regexpr(",",Data.Analyzed)!=-1, Data.Analyzed:='Joint Analysis']
  } else out.dt[, Data.Analyzed:=this.reference]
    
  out.dt[, Compound:=fitdata$Compound[1]]

  # If any of the parameters were not optimized or if the the model does not fit well:
#  if (any(ln.means==as.vector(opt.params[names(ln.means)])) | any(sigmas>1))
  if (any(ln.means==as.vector(orig.params[names(ln.means)])) | any(sigmas>MAXSIGMA))
  {
    out.dt[,LogLikelihood:=0]
    out.dt[,AIC:=Inf]
    cat(paste("Some parameters were not optimized or sigma >",MAXSIGMA,". Returning AIC=INF.\n"))
  } else if (model=='1compartment')
  {
    # Check for bad one-compartment model fits:
    if (ln.means["Vdist"]==upper["Vdist"])
    {
      out.dt[,LogLikelihood:=0]
      out.dt[,AIC:=Inf]
      cat("Vdist equal to upper bound on optimizer. Returning AIC=INF.\n")
    } else {
      out.dt[,LogLikelihood:=-all.data.fit$value]
      out.dt[,AIC:=2*length(opt.params)+2*all.data.fit$value]
    }
  } else if (model=='2compartment')
  {
    # Check for bad two-compartment model fits:
    if (ln.means["Fbetaofalpha"]==upper["Fbetaofalpha"]|ln.means["Fbetaofalpha"]==lower["Fbetaofalpha"]|ln.means["V1"]==upper["V1"])
    {
      out.dt[,LogLikelihood:=0]
      out.dt[,AIC:=Inf]
      cat("Problem with Fbetaofalpha or V1 equaling optimizer bound. Returning AIC=INF.\n")
    } else {
      out.dt[,LogLikelihood:=-all.data.fit$value]
      out.dt[,AIC:=2*length(opt.params)+2*all.data.fit$value]
    }
  }
  
  if (!is.null(this.reference))
  {
    out.dt[,Reference:=NULL]
  }

  return(out.dt)
}

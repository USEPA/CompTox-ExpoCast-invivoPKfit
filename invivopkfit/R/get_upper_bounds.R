get_upper_bounds <- function(fitdata,
                             this.dtxsid,
                             modelfun,
                             model,
                             this.reference = NULL,
                             suppress.messages = FALSE,
                             sig.figs = 5,
                             UPPERBOUNDARY = 1e4){
  #specify upper bounds of params to optimize (on a log scale!!)
  upper[] <- log(UPPERBOUNDARY)
  TRYSIGMA <- stats::median(fitdata$Value, na.rm = TRUE)
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

}

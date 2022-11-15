analyze_subset <- function(fitdata,
                           this.dtxsid,
                           modelfun,
                           model,
                           this.reference = NULL,
                           suppress.messages = FALSE,
                           sig.figs = 5){



  #log-transform the model parameters
  these.params <- lapply(these.params, log)

  #Add the per-study standard deviation to the set of params to optimize
  opt.params <- c(opt.params,
                  these.params[regexpr("sigma2",
                                       names(these.params)) != -1])

  MAXSIGMA <- 100 * stats::median(fitdata$Value, na.rm = T)
}

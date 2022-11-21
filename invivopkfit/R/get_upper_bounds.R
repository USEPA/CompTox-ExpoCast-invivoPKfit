#' Get upper bounds for estimating model parameters
#'
#' For a set of model parameters, get the upper bounds for the optimizer (if
#' parameter is to be estimated).
#'
#' @param fitdata A data.frame: the concentration-time-dose data to be used for
#'   fitting.
#' @param par_DF Optional: A data.frame as produced by [get_opt_params], with a
#'   character variable [param_name] containing parameter names, and a logical
#'   variable [optimize_param] containing TRUE if parameter is to be fitted, and
#'   FALSE if parameter is to be held constant. Any other variables will be
#'   ignored. Default NULL, in which case it will be determined by calling
#'   [get_opt_params].
#' @param model The name of the model whose parameters are to be estimated.
#'   Currently only "flat", "1compartment", or "2compartment" is supported.
#'   Ignored if `par_DF` is provided.
#' @param pool_sigma Logical: Whether to pool all data (estimate only one error
#'   standard deviation) or not (estimate separate error standard deviations for
#'   each reference). Default FALSE to estimate separate error SDs for each
#'   reference. (If `fitdata` only includes one reference, `pool_sigma` will
#'   have no effect.) Ignored if `par_DF` is provided.
#' @return A data.frame: `parDF` with additional variables `upper_bound`
#'   (numeric, containing the upper bound for each parameter) and
#'   `upper_bound_msg` (character, containing a brief message explaining how the
#'   upper-bound value was calculated).
#' @author Caroline Ring, John Wambaugh, Mitchell Teague
get_upper_bounds <- function(fitdata,
                             par_DF = NULL,
                             model,
                             pool_sigma = FALSE,
                             suppress.messages = FALSE){

  if(is.null(par_DF)){
  par_DF <- get_opt_params(model = model,
                           fitdata = fitdata,
                           pool_sigma = pool_sigma,
                           param_names = par_DF$param_name,
                           suppress.messages = suppress.messages)
  }
  rownames(par_DF) <- par_DF$param_name


  #Set a plausible upper bound for sigmas:
  #by default, this is 100 times the median of the concentrations.
  MAXSIGMA <- 100*stats::median(fitdata$Value, na.rm = TRUE)
  #for comparison, check the sd of the concentrations.
  sigma_test <- sd(fitdata$Value, na.rm = TRUE)

  #Catch cases where the sd of the data is much larger than the median.
  if (sigma_test > MAXSIGMA) {
    sigma_msg <- paste("BIGSD: Conc. std. dev., ",
                 format(sigma_test, digits = digits, scientific = scientific),
                 "is much higher than than usual upper bound = 100 * conc. median = ",
                 format(MAXSIGMA, digits = digits, scientific = scientific),
                 ". Using upper bound = 2 * SD = ",
                 format(2*sigma_test, digits = digits, scientific = scientific),
                 ".", sep="")
    if(!suppress.messages){
    message(sigma_msg)
    }
    big.sd <- TRUE
    MAXSIGMA <- 2*sigma_test
  } else {
    big.sd <- FALSE
    sigma_msg <- paste0("Upper bound = 100 * conc. median = ",
    format(MAXSIGMA, digits = digits, scientific = scientific))
    }

#Assign sigma upper bounds and upper-bound messages
  par_DF[grepl(x = par_DF$param_name,
               pattern = "sigma"),
         "upper_bound"] <- MAXSIGMA
  par_DF[grepl(x = par_DF$param_name,
               pattern = "sigma"),
         "upper_bound_msg"] <- sigma_msg

  #kelim
  par_DF[grepl(x = par_DF$param_name,
               pattern = "kelim"),
         c("upper_bound",
           "upper_bound_msg")] <- list(1e4,
                                       "Default")

  #Vdist or V1
  par_DF[grepl(x = par_DF$param_name,
               pattern = "Vdist|V1"),
         c("upper_bound",
           "upper_bound_msg")] <- list(max(fitdata$Dose) / min(fitdata$LOQ),
                                       "Max dose/min LOQ")

  #Ralphatokelim
  par_DF[grepl(x = par_DF$param_name,
               pattern = "Ralphatokelim"),
         c("upper_bound",
           "upper_bound_msg")] <- list(1000,
                                       "Default")

  #Fbetaofoalpha
  par_DF[grepl(x = par_DF$param_name,
               pattern = "Fbetaofalpha"),
         c("upper_bound",
           "upper_bound_msg")] <- list(0.75,
                                       "Default")

  #Fgutabs
  par_DF[grepl(x = par_DF$param_name,
               pattern = "Fgutabs"),
         c("upper_bound",
           "upper_bound_msg")] <- list(1,
                                       "Default")

  #kgutabs
  par_DF[grepl(x = par_DF$param_name,
               pattern = "kgutabs"),
         c("upper_bound",
           "upper_bound_msg")] <- list(1000,
                                       "Default")

  #For anything not to be optimized, set its bounds to NA
  par_DF[!(par_DF$optimize_param %in% TRUE),
         c("upper_bound",
           "upper_bound_msg")] <- list(NA_real_,
                                       "optimize_param is not TRUE")

  #For anything not to be used, set its bounds to NA
  par_DF[!(par_DF$use_param %in% TRUE),
         c("upper_bound",
           "upper_bound_msg")] <- list(NA_real_,
                                       "use_param is not TRUE")

  return(par_DF)

}

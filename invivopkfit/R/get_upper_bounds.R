#' Get upper bounds for fitting model parameters
#'
#' For a set of model parameters, get the upper bounds for the optimizer (if
#' parameter is to be fitted).
#'
#' @param par_DF Optional: A data.frame as produced by [get_opt_params], with a
#'   character variable [param_name] containing parameter names, and a logical
#'   variable [optimize_param] containing TRUE if parameter is to be fitted, and
#'   FALSE if parameter is to be held constant. Any other variables will be
#'   ignored. Default NULL, in which case it will be determined by calling
#'   [get_opt_params].
#' @param model The name of the model whose parameters are to be estimated.
#' @param fitdata A data.frame: the concentration-time-dose data to be used for
#'   fitting.
#' @return A data.frame: `parDF` with additional variables `upper_value`
#'   (numeric, containing the upper bound for each parameter) and
#'   `upper_value_msg` (character, containing a brief message explaining how the
#'   upper-bound value was calculated).
#' @author Caroline Ring, John Wambaugh, Mitchell Teague
get_upper_bounds <- function(fitdata,
                             model,
                             par_DF = NULL,
                             digits = 3,
                             scientific = -2,
                             UPPERBOUNDARY = 1e4){

  if(is.null(par_DF)){
  par_DF <- get_opt_params(model = model,
                           fitdata = fitdata,
                           param_names = par_DF$param_name)
  }
  rownames(par_DF) <- par_DF$param_name


  #Set a plausible upper bound for sigmas:
  #by default, this is the median of the concentrations.
  MAXSIGMA <- stats::median(fitdata$Value, na.rm = TRUE)
  #for comparison, check the sd of the concentrations.
  sigma_test <- sd(fitdata$Value, na.rm = TRUE)

  #Catch cases where the sd of the data is much larger than the median.
  if (sigma_test > MAXSIGMA) {
    sigma_msg <- paste("BIGSD: Conc. std. dev., ",
                 format(sigma_test, digits = digits, scientific = scientific),
                 "is much higher than than usual upper bound = conc. median = ",
                 format(MAXSIGMA, digits = digits, scientific = scientific),
                 ". Using upper bound = 2 * SD = ",
                 format(2*sigma_test, digits = digits, scientific = scientific),
                 ".", sep="")
    warning(msg)
    big.sd <- TRUE
    MAXSIGMA <- 2*sigma_test
  } else {
    big.sd <- FALSE
    sigma_msg <- paste0("Upper bound = conc. median = ",
    format(MAXSIGMA, digits = digits, scientific = scientific))
    }

#Assign sigma upper bounds and upper-bound messages
  par_DF[grepl(x = par_DF$param_name,
               pattern = "sigma"),
         "upper_value"] <- MAXSIGMA
  par_DF[grepl(x = par_DF$param_name,
               pattern = "sigma"),
         "upper_value_msg"] <- sigma_msg

  #kelim
  par_DF[grepl(x = par_DF$param_name,
               pattern = "kelim"),
         c("upper_value",
           "upper_value_msg")] <- list(1e4,
                                       "Default")

  #Vdist or V1
  par_DF[grepl(x = par_DF$param_name,
               pattern = "Vdist|V1"),
         c("upper_value",
           "upper_value_msg")] <- list(max(fitdata$Dose) / min(fitdata$LOQ),
                                       "Max dose/min LOQ")

  #Ralphatokelim
  par_DF[grepl(x = par_DF$param_name,
               pattern = "Ralphatokelim"),
         c("upper_value",
           "upper_value_msg")] <- list(1000,
                                       "Default")

  #Fbetaofoalpha
  par_DF[grepl(x = par_DF$param_name,
               pattern = "Fbetaofalpha"),
         c("upper_value",
           "upper_value_msg")] <- list(0.75,
                                       "Default")

  #Fgutabs

  par_DF[grepl(x = par_DF$param_name,
               pattern = "Fgutabs"),
         c("upper_value",
           "upper_value_msg")] <- list(1000,
                                       "Default")

  #kgutabs
  par_DF[grepl(x = par_DF$param_name,
               pattern = "kgutabs"),
         c("upper_value",
           "upper_value_msg")] <- list(1,
                                       "Default")

  #For anything not to be optimized, set its bounds to NA
  par_DF[!(optimize_param %in% TRUE),
         c("upper_value",
           "upper_value_msg")] <- list(NA_real_,
                                       "optimize_param is not TRUE")

  return(par_DF)

}

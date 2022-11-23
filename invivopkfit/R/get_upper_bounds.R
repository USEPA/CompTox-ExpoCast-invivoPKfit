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
#' @param upper_default A `data.frame` with two variables: `param_name`, giving
#'   the names of all parameters for all models, plus an additional parameter
#'   "sigma" (representing the upper bound applied to all residual standard
#'   deviations); `upper_bound`, giving the default upper-bound values for each
#'   parameter; and `upper_bound_msg`, giving a message about the default upper
#'   bound values.
#' @param sigma_upper_from_data Logical: if TRUE, override the value for "sigma"
#'   in `upper_default` and instead use either 100 times the median
#'   concentration, or 2 times the overall standard deviation of concentrations,
#'   whichever is greater. If FALSE, use the upper bound for "sigma" in
#'   `upper_default.` Default value is TRUE.
#' @param Vdist_upper_from_data Logical: if TRUE, override the value for "Vdist"
#'   or "V1" in `upper_default`, and instead use `Vdist_factor` times the
#'   maximum dose in the data, divided by the minimum concentration or LOQ in
#'   the data. If FALSE, use the upper bound for "Vdist" or "V1" in
#'   `upper_default`. Default value is TRUE.
#' @param Vdist_factor Numeric: If `Vdist_upper_from_data` is TRUE, then
#'   `Vdist_factor` scales the ratio between the maximum dose and minimum
#'   concentration to arrive at a theoretical upper bound for the volume of
#'   distribution ("vdist" or "V1"). Default value is 1. If
#'   `Vdist_upper_from_data` is FALSE, then `Vdist_factor` has no effect.
#' @return A data.frame: `parDF` with additional variables `upper_bound`
#'   (numeric, containing the upper bound for each parameter) and
#'   `upper_bound_msg` (character, containing a brief message explaining how the
#'   upper-bound value was calculated).
#' @author Caroline Ring, John Wambaugh, Mitchell Teague
get_upper_bounds <- function(fitdata,
                             par_DF = NULL,
                             model,
                             pool_sigma = FALSE,
                             upper_default = data.frame(
                               param_name = c("A",
                                              "kelim",
                                              "Vdist",
                                              "kgutabs",
                                              "Fgutabs",
                                              "V1",
                                              "k12",
                                              "k21",
                                              "Fgutabs_Vdist",
                                              "Fgutabs_V1",
                                              "sigma"),
                               upper_bound = c(1e8, #A
                                               1e4, #kelim
                                               1e8, #Vdist
                                               1e4, #kgutabs
                                               1, #Fgutabs
                                               1e8, #V1
                                               1e4, #k12
                                               1e4, #k21
                                               1e8, #Fgutabs_Vdist
                                               1e8, #Fgutabs_V1
                                               1e8), #sigma
                               upper_bound_msg = "Default"
                             ),
                             sigma_upper_from_data = TRUE,
                             Vdist_upper_from_data = TRUE,
                             Vdist_factor = 2,
                             suppress.messages = FALSE,
                             digits = 5,
                             scientific = -3){

  if(is.null(par_DF)){
  par_DF <- get_opt_params(model = model,
                           fitdata = fitdata,
                           pool_sigma = pool_sigma,
                           param_names = par_DF$param_name,
                           suppress.messages = suppress.messages)
  }
  rownames(par_DF) <- par_DF$param_name

  #Use the default upper bounds
  #this gets everything except sigma, which will not be in the model params
  par_DF <- merge(par_DF,
                  upper_default,
                  by = "param_name",
                  all.x = TRUE,
                  all.y = FALSE)

  if(sigma_upper_from_data %in% TRUE){
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

  par_DF[grepl(x = par_DF$param_name,
               pattern = "sigma"),
         c("upper_bound", "upper_bound_msg")] <- list(MAXSIGMA,
                                                      sigma_msg)
  }else{
    #use defaults
    par_DF[grepl(x = par_DF$param_name,
                 pattern = "sigma"),
           "upper_bound"] <- upper_default[upper_default$param_name %in% "sigma",
                                           "upper_bound"]
    par_DF[grepl(x = par_DF$param_name,
                 pattern = "sigma"),
           "upper_bound_msg"] <- upper_default[upper_default$param_name %in% "sigma",
                                               "upper_bound_msg"]
  }

  if(Vdist_upper_from_data %in% TRUE){
  #Overwrite Vdist/V1 max with a value taken from the data
  par_DF[grepl(x = par_DF$param_name,
               pattern = "Vdist|V1"),
         c("upper_bound",
           "upper_bound_msg")] <- list(Vdist_factor *
                                         max(fitdata$Dose) /
                                         min(
                                           pmin(fitdata$LOQ,
                                                                  fitdata$Value,
                                                                  na.rm = TRUE)
                                           ),
                                       paste(Vdist_factor,
                                             "* max dose/[min conc or LOQ]")
           )
  }

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

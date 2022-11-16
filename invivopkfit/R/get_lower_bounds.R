#' Get lower bounds for fitting model parameters
#'
#' For a set of model parameters, get the lower bounds for the optimizer (if
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
#' @return A data.frame: `parDF` with additional variables `lower_value`
#'   (numeric, containing the lower bound for each parameter) and
#'   `lower_value_msg` (character, containing a brief message explaining how the
#'   lower-bound value was calculated).
#' @author Caroline Ring, John Wambaugh, Mitchell Teague
#'
get_lower_bounds <- function(par_DF = NULL,
                             model,
                             fitdata){
  if(is.null(par_DF)){
    par_DF <- get_opt_params(model = model,
                             fitdata = fitdata,
                             param_names = par_DF$param_name)
  }
  rownames(par_DF) <- par_DF$param_name

  #defaults
  par_DF[, c("lower_value",
             "lower_value_msg")] <- list(1e-8, "Default")

  #Override default for params where it doesn't make sense

  #Vdist/V1
  par_DF[grepl(x = par_DF$param_names,
               pattern = "Vdist|V1"),
         c("lower_value",
           "lower_value_msg")] <- list(0.01,
                                       "Default")

  #Ralphatokelim
  par_DF[grepl(x = par_DF$param_names,
               pattern = "Ralphatokelim"),
         c("lower_value",
           "lower_value_msg")] <- list(0.0,
                                       "Default")

  #Fgutabs
  par_DF[grepl(x = par_DF$param_names,
               pattern = "Fgutabs"),
         c("lower_value",
           "lower_value_msg")] <- list(0.05,
                                       "Default")

  #sigma
  par_DF[grepl(x = par_DF$param_names,
               pattern = "sigma"),
         c("lower_value",
           "lower_value_msg")] <- list(1e-5,
                                       "Default")

  #For anything not to be optimized, set its bounds to NA
  par_DF[!(optimize_param %in% TRUE),
         c("lower_value",
           "lower_value_msg")] <- list(NA_real_,
                                       "optimize_param is not TRUE")

  return(par_DF)
}

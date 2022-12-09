#' Get upper bounds for estimating model parameters
#'
#' For a set of model parameters, get the upper bounds for the optimizer (if
#' parameter is to be estimated).
#'
#' The default `upper_default` `data.frame` is shown below in table format:
#'
#' | param_name     | upper_bound | upper_bound_msg |
#' | ---------------| ----------- | --------------- |
#' | A              | Inf           | Default         |
#' | kelim          | Inf        | Default         |
#' | Vdist          | Inf        | Default         |
#' | kgutabs        | Inf        | Default         |
#' | Fgutabs        | 1         | Default         |
#' | V1             | Inf           | Default         |
#' | k12            | Inf         | Default         |
#' | k21            | Inf         | Default         |
#' | Fgutabs_Vdist  | Inf    | Default         |
#' | Fgutabs_V1     | Inf    | Default         |
#' | sigma          | Inf         | Default         |
#'
#' @param fitdata A data.frame: the concentration-time-dose data to be used for
#'   fitting.
#' @param par_DF Optional: A data.frame as produced by [get_opt_params()], with a
#'   character variable `param_name` containing parameter names, and a logical
#'   variable `optimize_param` containing TRUE if parameter is to be fitted, and
#'   FALSE if parameter is to be held constant. Any other variables will be
#'   ignored. Default NULL, in which case it will be determined by calling
#'   [get_opt_params()].
#' @param model The name of the model whose parameters are to be estimated.
#'   Currently only "flat", "1compartment", or "2compartment" is supported.
#'   Ignored if `par_DF` is provided.
#' @param pool_sigma Logical: Whether to pool all data (estimate only one error
#'   standard deviation) or not (estimate separate error standard deviations for
#'   each reference). Default FALSE to estimate separate error SDs for each
#'   reference. (If `fitdata` only includes one reference, `pool_sigma` will
#'   have no effect.) Ignored if `par_DF` is provided.
#' @param upper_default A `data.frame` with three variables: `param_name`,
#'   giving the names of parameters; `upper_bound`, giving the default
#'   upper-bound values for each parameter; and `upper_bound_msg`, giving a
#'   message about the default upper bound values. See Details for default
#'   value.
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
                               upper_bound = c(Inf, #A
                                               Inf, #kelim
                                               Inf, #Vdist
                                               Inf, #kgutabs
                                               1, #Fgutabs
                                               Inf, #V1
                                               Inf, #k12
                                               Inf, #k21
                                               Inf, #Fgutabs_Vdist
                                               Inf, #Fgutabs_V1
                                               Inf), #sigma
                               upper_bound_msg = "Default"
                             ),
                             sigma_upper_from_data = TRUE,
                             suppress.messages = FALSE){

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
  #It really shouldn't ever be that much worse than the null model...
    Astart <- exp(mean(log(fitdata$Value/fitdata$Dose), na.rm = TRUE))
    log_resid_flat <-  log(Astart * fitdata$Dose) - log(fitdata$Value)
    log_resid_flat[!is.finite(log_resid)] <- NA_real_
    sd_flat <- sd(log_resid_flat, na.rm = TRUE)

#...but give it a factor of 10 just to be sure
  par_DF[grepl(x = par_DF$param_name,
               pattern = "sigma"),
         c("upper_bound", "upper_bound_msg")] <- list(10*sd_flat,
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

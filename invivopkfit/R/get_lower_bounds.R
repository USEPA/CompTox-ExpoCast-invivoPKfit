#' Get lower bounds for estimating model parameters
#'
#' For a set of model parameters, get the lower bounds for the optimizer.
#'
#' The default `lower_default` `data.frame` is shown below in table format:
#'
#' | param_name     | lower_bound | lower_bound_msg |
#' | ---------------| ----------- | --------------- |
#' | kelim          | 1e-4        | Default         |
#' | Vdist          | 1e-4        | Default         |
#' | kgutabs        | 1e-4        | Default         |
#' | Fgutabs        | 1e-4         | Default         |
#' | V1             | 1e-4           | Default         |
#' | k12            | 1e-4         | Default         |
#' | k21            | 1e-4         | Default         |
#' | Fgutabs_Vdist  | 1e-4    | Default         |
#' | Fgutabs_V1     | 1e-4    | Default         |
#' | sigma          | 1e-8         | Default         |
#'
#'
#' Any parameters which will not be estimated from data (based on either the
#' variable `optimize_param` in `par_DF` if `par_DF` is provided, or the output
#' of [get_opt_params()] if `par_DF` is not provided) are assigned a lower bound
#' of `NA_real_`.
#'
#' @param fitdata A `data.frame`: the concentration-time-dose data to be used for
#'   fitting, for example as produced by [preprocess_data()].
#' @param par_DF Optional: A `data.frame` as produced by [get_opt_params()], with a
#'   character variable `param_name` containing parameter names, and a logical
#'   variable `optimize_param` containing `TRUE` if parameter is to be fitted, and
#'   `FALSE` if parameter is to be held constant. Any other variables will be
#'   ignored. Default `NULL`, in which case it will be determined by calling
#'   [get_opt_params()].
#' @param model The name of the model whose parameters are to be estimated.
#'   Currently only "flat", "1compartment", or "2compartment" is supported.
#'   Ignored if `par_DF` is provided.
#' @param pool_sigma Logical: Whether to pool all data (estimate only one error
#'   standard deviation) or not (estimate separate error standard deviations for
#'   each study). Default `FALSE` to estimate separate error SDs for each
#'   study. (If `fitdata` only includes one study, `pool_sigma` will
#'   have no effect.) Ignored if `par_DF` is provided.
#' @param lower_default  A `data.frame` with three variables: `param_name`,
#'   giving the names of parameters; `lower_bound`, giving the default
#'   lower-bound values for each parameter; and `lower_bound_msg`, giving a
#'   message about the default lower bound values. See Details for default
#'   value.
#' @return A `data.frame`: the same as `parDF` with additional variables `lower_bound`
#'   (numeric, containing the lower bound for each parameter) and
#'   `lower_bound_msg` (character, containing a brief message explaining how the
#'   lower-bound value was calculated).
#' @author Caroline Ring, John Wambaugh, Mitchell Teague
#'
get_lower_bounds <- function(fitdata,
                             par_DF = NULL,
                             model,
                             pool_sigma = FALSE,
                             lower_default = data.frame(
                               param_name = c("kelim",
                                              "Vdist",
                                              "kgutabs",
                                              "Fgutabs",
                                              "V1",
                                              "k12",
                                              "k21",
                                              "Fgutabs_Vdist",
                                              "Fgutabs_V1",
                                              "Rblood2plasma",
                                              "sigma"),
                               lower_bound = c(1e-4, #kelim
                                               1e-4, #Vdist
                                               1e-4, #kgutabs
                                               1e-4, #Fgutabs
                                               1e-4, #V1
                                               1e-4, #k12
                                               1e-4, #k21
                                               1e-4, #Fgutabs_Vdist
                                               1e-4, #Fgutabs_V1
                                               1e-4, #Rblood2plasma
                                               1e-8 #sigma
                                               ),
                               lower_bound_msg = "Default"
                             ),
                             suppress.messages = FALSE
                             ){
  if(is.null(par_DF)){
    par_DF <- get_opt_params(model = model,
                             fitdata = fitdata,
                             pool_sigma = pool_sigma,
                             param_names = par_DF$param_name,
                             suppress.messages = suppress.messages)
  }
  rownames(par_DF) <- par_DF$param_name

  #Use the default lower bounds
  #this gets everything except sigma, which will not be in the model params
  par_DF <- merge(par_DF,
                  lower_default,
                  by = "param_name",
                  all.x = TRUE,
                  all.y = FALSE)


  #sigma
  par_DF[grepl(x = par_DF$param_name,
               pattern = "sigma"),
         c("lower_bound",
           "lower_bound_msg")] <- lower_default[
             lower_default$param_name %in% "sigma",
             c("lower_bound",
               "lower_bound_msg")
           ]
#set rownames to param names
  rownames(par_DF) <- par_DF$param_name

  #For anything not to be optimized, set its bounds to NA
  #because the bounds will not be used
  par_DF[!(par_DF$optimize_param %in% TRUE),
         c("lower_bound",
           "lower_bound_msg")] <- list(NA_real_,
                                       "optimize_param is not TRUE")

  #For anything not to be used, set its bounds to NA
  par_DF[!(par_DF$use_param %in% TRUE),
         c("lower_bound",
           "lower_bound_msg")] <- list(NA_real_,
                                       "use_param is not TRUE")

  return(par_DF)
}

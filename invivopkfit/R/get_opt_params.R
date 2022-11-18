#' Determine which parameters to fit (optimize)
#'
#' Determine whether to fit (optimize) each model parameter from the data or
#' hold it constant
#'
#' All model parameters will be fitted (optimized) with the following
#' exceptions:
#'
#' For 1-compartment and 2-compartment models:
#'
#' If no oral dosing data are available, then the parameter "kgutabs" will not
#' be estimated from the data. Instead it will be held constant at its
#' "starting" value while other parameters are estimated.
#'
#' If data do not include both oral dosing and IV dosing, then the parameter
#' "Fgutabs" will not be estimated from the data. Instead it will be held
#' constant at its "starting" value while other parameters are estimated.
#'
#' If non-NULL `param_names` do not match the result of
#' `get_model_parameters(param_names)`, this function will throw a warning, but
#' will otherwise proceed using `param_names`.
#'
#' @param model Character: the name of the model to parameterize.
#' @param fitdata Data.frame: The set of concentration-dose-time data to be used
#'   to fit the model parameters. Requires a variable named "Route" which
#'   contains the dosing route, either "po" (oral dosing) or "iv" (intravenous
#'   dosing).
#' @param param_names Optional: A character vector of parameter names needed by
#'   the model. Default NULL to automatically determine these from `model` by
#'   calling [get_model_paramnames].
#'
#' @return A named logical vector whose names are `param_names`, indicating
#'   whether to fit each model parameter (TRUE) or hold it constant (FALSE).
#'
#'

get_opt_params <- function(model,
                           fitdata,
                           param_names = NULL,
                           sigma_ref = TRUE,
                           suppress.messages = FALSE){

  model_param_names <- get_model_paramnames(model = model)
  if(is.null(param_names)){
  param_names <- model_param_names
  }else{
    if(!all(model_param_names %in% param_names) |
       !all(param_names %in% model_param_names)){
      warning(paste0("get_opt_params(): param_names do not match the result of ",
      "`get_model_paramnames(model).\n",
      "param_names = ",
      paste(param_names, collapse = ", "),
      "\n",
      "get_model_paramnames(model) = ",
      paste(model_param_names, collapse = ", "),
      "\n",
      "Proceeding using param_names and ignoring model."
      ))
    }
  }

  opt_params <- rep(TRUE, length(param_names))
  names(opt_params) <- param_names

  if("kgutabs" %in% param_names &
     !("po" %in% fitdata$Route)){
    #if no oral data, can't fit kgutabs
    opt_params["kgutabs"] <- FALSE
  }

  if("Fgutabs" %in% param_names &
     !(all(c("po", "iv") %in% fitdata$Route))){
    #if we don't have both IV and oral data,
    #can't fit Fgutabs
    opt_params["Fgutabs"] <- FALSE
  }

  par_DF <- data.frame(param_name = names(opt_params),
                       optimize_param = opt_params)

  #add hyperparameters: individual reference error SDs if sigma_ref = TRUE
  if(sigma_ref %in% TRUE){
 refs <- unique(fitdata$Reference)
 # Mark these as parameters to be optimized
 sigma_ref_DF <- data.frame(param_name = paste("sigma",
                                               "ref",
                                               refs,
                                               sep = "_"),
                            optimize_param = TRUE)
  }else{ #just add a general "sigma" parameter, not one for each individual reference
    sigma_ref_DF <- data.frame(param_name = "sigma",
                               optimize_param = TRUE)
  }

  par_DF <- rbind(par_DF, sigma_ref_DF)

  return(par_DF)
}

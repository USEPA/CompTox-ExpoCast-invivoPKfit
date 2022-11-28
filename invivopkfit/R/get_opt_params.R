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
#' @param pool_sigma Logical: Whether to pool all data (estimate only one error
#'   standard deviation) or not (estimate separate error standard deviations for
#'   each reference). Default FALSE to estimate separate error SDs for each
#'   reference. (If `fitdata` only includes one reference, `pool_sigma` will have
#'   no effect.)
#'
#' @return A `data.frame` with the following variables:
#' - `param_name` Character: the name of each model parameter.
#' - `optimize_param` Logical: Whether the parameter is to be
#' optimized/estimated from the data (TRUE), or not (FALSE). For example, if
#' `model` is '1compartment' and `fitdata` does not contain both oral and IV
#' data, then the parameter "Fgutabs" cannot be estimated from the data, so
#' `optimize_param` for "Fgutabs" will be FALSE.
#' - `use_param` Logical: Whether the parameter gets used in the model at all
#' (TRUE) or not (FALSE). For example, if `model` is '1compartment' and
#' `fitdata` contains no oral data at all, then for parameter "Fgutabs",
#' `use_param` will be FALSE, because the 1-compartment model function for
#' IV-only data will not use the value for "Fgutabs".
#'
#'  By default, `use_param` and `opt_param` will be the same -- for each
#'  parameter, both `TRUE` or both `FALSE`. However, both variables are retained
#'  to allow for the possibility of holding one or more parameters constant and
#'  optimizing the rest (which would be accomplished by setting
#' `optimize_param  = FALSE` but keeping `use_param = TRUE` for the parameters
#' to be held constant).
#'
#'

get_opt_params <- function(model,
                           fitdata,
                           param_names = NULL,
                           pool_sigma = FALSE,
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

  use_params <- rep(TRUE, length(param_names))
  names(use_params) <- param_names

  if("kgutabs" %in% param_names &
     !("po" %in% fitdata$Route)){
    #if no oral data, can't fit kgutabs,
    #and it won't be used.
    opt_params["kgutabs"] <- FALSE
    use_params["kgutabs"] <- FALSE
  }

  if("Fgutabs" %in% param_names &
     !("po" %in% fitdata$Route)){
    #if no oral data, can't fit Fgutabs,
    #and it won't be used.
    opt_params["Fgutabs"] <- FALSE
    use_params["Fgutabs"] <- FALSE
  }

  if("Fgutabs" %in% param_names &
     ("po" %in% fitdata$Route)){
    #if there *is* oral data
    #if oral data...
    if("iv" %in% fitdata$Route){
      #if both oral and IV data,
      #then fit Fgutabs and Vdist separately,
      #so turn off Fgutabs_Vdist
      if("Fgutabs_Vdist" %in% param_names){
      opt_params["Fgutabs_Vdist"] <- FALSE
      use_params["Fgutabs_Vdist"] <- FALSE
      }
      if("Fgutabs_V1" %in% param_names){
        opt_params["Fgutabs_V1"] <- FALSE
        use_params["Fgutabs_V1"] <- FALSE
      }
    }else{
      #if oral but no IV data,
      #then turn off Fgutabs and Vdist,
      #and just fit Fgutabs_Vdist
      opt_params["Fgutabs"] <- FALSE
      use_params["Fgutabs"] <- FALSE
      if("Vdist" %in% param_names){
        opt_params["Vdist"] <- FALSE
        use_params["Vdist"] <- FALSE
      }
      if("V1" %in% param_names){
        opt_params["V1"] <- FALSE
        use_params["V1"] <- FALSE
      }
    }
  }

  par_DF <- data.frame(param_name = names(opt_params),
                       optimize_param = opt_params,
                       use_param = use_params)

  #add hyperparameters: individual reference error SDs
  refs <- unique(fitdata$Reference)

  #if more than one reference and user has not specified to pool data:
  if(pool_sigma %in% FALSE &
     length(refs) > 1){
    # Add individual reference error SDs, named as sigma_ref_[Reference ID]
    # Mark these as parameters to be optimized
 sigma_ref_DF <- data.frame(param_name = paste("sigma",
                                               "ref",
                                               refs,
                                               sep = "_"),
                            optimize_param = TRUE,
                            use_param = TRUE)
  }else{ #just add a single pooled "sigma" parameter
    sigma_ref_DF <- data.frame(param_name = "sigma",
                               optimize_param = TRUE,
                               use_param = TRUE)
  }

  par_DF <- rbind(par_DF, sigma_ref_DF)

  return(par_DF)
}

#' Determine which parameters to optimize from data
#'
#' Determine whether to optimize each model parameter from the data
#'
#' Typically this function is not called directly by the user. Rather, it is
#' called by [analyze_subset()] which in turn is called by the main fitting
#' function, [fit_all()].
#'
#' It may not be possible to estimate the full set of model parameters from a
#' given dataset, depending on which routes of dose administration are
#' represented in the dataset (oral and/or intravenous dosing). Some parameters
#' can be identified only if oral dosing data are available. Some can be
#' identified only if both oral and intravenous dosing data are available. This
#' function determines which parameters will and will not be optimized, based on
#' whether oral and/or intravenous dosing data are available.
#'
#' # 1-compartment model
#'
#' The full set of model parameters for the 1-compartment model includes
#' `Vdist`, `kelim`, `kgutabs`, and `Fgutabs`.
#'
#' ## IV data, no oral data
#'
#' If IV dosing data are available, but no oral dosing data are available, then
#' only the parameters `Vdist` and `kelim` will be estimated from the data. The
#' parameters `kgutabs` and `Fgutabs` cannot be estimated from IV data alone.
#'
#' ## Oral data, no IV data
#'
#' If oral dosing data are available, but no IV dosing data are available, then
#' the parameters `kelim` and `kgutabs` can be estimated from the data. However,
#' the parameters `Fgutabs` and `Vdist` cannot be identified separately. From
#' oral data alone, only the ratio `Fgutabs/Vdist` can be identified. This ratio
#' is represented by a single parameter named `Fgutabs_Vdist`. `Fgutabs` and
#' `Vdist` will not be optimized, but `Fgutabs_Vdist` will be optimized, along
#' with `kelim` and `kgutabs`.
#'
#' ## Oral data and IV data
#'
#' If both oral and IV dosing data are available, then `Vdist`, `kelim`,
#' `kgutabs`, and `Fgutabs` will all be estimated from the data.
#'
#' # 2-compartment model
#'
#' The full set of model parameters for the 1-compartment model includes `V1`,
#' `kelim`, `k12`, `k21`, `kgutabs`, and `Fgutabs`.
#'
#' ## IV data, no oral data
#'
#' If IV dosing data are available, but no oral dosing data are available, then
#' only the parameters `V1`, `kelim`, `k12`, and `k21` will be estimated from
#' the data. The parameters `kgutabs` and `Fgutabs` cannot be estimated from IV
#' data alone.
#'
#' ## Oral data, no IV data
#'
#' If oral dosing data are available, but no IV dosing data are available, then
#' the parameters `kelim`, `k12`, `k21`, and `kgutabs` will be estimated from
#' the data. However, the parameters `Fgutabs` and `V1` cannot be identified
#' separately. From oral data alone, only the ratio `Fgutabs/V1` can be
#' identified. This ratio is represented by a single parameter named
#' `Fgutabs_V1`. `Fgutabs` and `V1` will not be optimized, but `Fgutabs_V1` will
#' be optimized, along with `kelim`, `k12`, `k21`, and `kgutabs`.
#'
#' ## Oral data and IV data
#'
#' If both oral and IV dosing data are available, then `V1`, `kelim`, `k12`,
#' `k21`, `kgutabs`, and `Fgutabs` will all be estimated from the data.
#'
#' # Further note on parameter identifiability
#'
#' This function does *not* guarantee that all parameters to be estimated are in
#' fact identifiable from the available data. This function does not check that
#' a sufficient number of observations are available, nor that those
#' observations cover sufficient time, to identify all parameters.
#'
#' For example, in order to fully identify all parameters of the 2-compartment
#' model, concentration measurements must be made at time points that capture
#' all three phases: absorption, distribution, and elimination phases. If no
#' measurements were made late enough to capture the elimination phase, then
#' `kelim` may not be identifiable from the data. This function makes no attempt
#' to detect such situations.
#'
#' # If user-specified parameter names do not match user-specified model
#'
#' If user-specified non-NULL `param_names` do not match the result of
#' [get_model_parameters()] for the user-specified `model`, this function will
#' throw a warning, but will otherwise proceed using the user-specified
#' `param_names`.
#'
#' @param model Character: the name of the model to parameterize.
#' @param fitdata `data.frame`: The set of concentration-dose-time data to be
#'   used to fit the model parameters, for example as produced by
#'   [preprocess_data()]. Requires a variable named `Route` which contains the
#'   dosing route, either "po" (oral dosing) or "iv" (intravenous dosing); a
#'   variable named `Media` which contains the medium in which concentration
#'   was measured (either "blood" or "plasma"), and a variable named `Reference`
#'   which contains reference identifiers.
#' @param param_names Optional: A character vector of parameter names needed by
#'   the model. Default NULL to automatically determine these from `model` by
#'   calling [get_model_paramnames()].
#' @param pool_sigma Logical: Whether to pool all data (estimate only one error
#'   standard deviation) or not (estimate separate error standard deviations for
#'   each reference). Default FALSE to estimate separate error SDs for each
#'   reference. (If `fitdata` only includes one reference, `pool_sigma` will
#'   have no effect.)
#'
#' @return A `data.frame` with the following variables:
#' - `param_name` Character: the name of each model parameter.
#' - `optimize_param` Logical: Whether the parameter is to be
#'   optimized/estimated from the data (TRUE), or not (FALSE).
#' - `use_param` Logical: Whether the parameter gets used in the model at all
#'   (TRUE) or not (FALSE).
#'
#'   By default, `use_param` and `opt_param` will be the same -- for each
#'   parameter, both `TRUE` or both `FALSE`. However, both variables are
#'   retained to allow for the possibility of holding one or more parameters
#'   constant and optimizing the rest (which would be accomplished by setting
#'   `optimize_param  = FALSE` but keeping `use_param = TRUE` for the parameters
#'   to be held constant).
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

  if("Fgutabs_Vdist" %in% param_names &
     !("po" %in% fitdata$Route)){
    #if no oral data, can't fit Fgutabs_Vdist,
    #and it won't be used.
    opt_params["Fgutabs_Vdist"] <- FALSE
    use_params["Fgutabs_Vdist"] <- FALSE
  }

  if("Fgutabs_V1" %in% param_names &
     !("po" %in% fitdata$Route)){
    #if no oral data, can't fit Fgutabs_V1,
    #and it won't be used.
    opt_params["Fgutabs_V1"] <- FALSE
    use_params["Fgutabs_V1"] <- FALSE
  }

  if(("po" %in% fitdata$Route)){
    #if there *is* oral data
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
      if("Fgutabs" %in% param_names){
      opt_params["Fgutabs"] <- FALSE
      use_params["Fgutabs"] <- FALSE
      }
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

  #if both "blood" and "plasma" are not in data,
  #then Rblood2plasma will not be optimized
  if(!(all(c("blood", "plasma") %in% fitdata$Media))){
    opt_params["Rblood2plasma"] <- FALSE
  }

  #if no medium is "blood" then Rblood2plasma will not be used at all
  if(!("blood" %in% fitdata$Media)){
    use_params["Rblood2plasma"] <- FALSE
    #otherwise, if blood-only data, Rblood2plasma will be used but held constant at 1
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

  #set rownames equal to param_name
  rownames(par_DF) <- par_DF$param_name

  return(par_DF)
}

#'Create a new `pk` object
#'
#'[pk()] initializes a new `pk` object.
#'
#'
#'[pk()] is used to construct the initial `pk` object for analysis. It is almost
#'always followed by `+` to add steps to the workflow. For example, you could
#'use `pk(my_data) + stat_model(model = '1comp')` to set up for fitting a
#'1-compartment model.
#'
#'# The `pk` object
#'
#'A `pk` object consists of a set of concentration-dose-time data to be fitted,
#'and sets of instructions for steps in the analysis workflow:
#'
#'- settings for how to pre-process the data (harmonizing variable names, imputing missing data, calculating derived variables)
#'- scalings/transformations to be applied to the data
#'- settings for the numerical optimization algorithm to be used to fit any model
#'- optionally: which PK model(s) should be fitted to this dataset. (You do not have to fit any PK model if you don't want to; you can instead just set up the `pk` object with data, and do non-compartmental analysis on it.)
#'
#'The most basic `pk` object is a named list with the following elements:
#'
#' - `data_orig`: The original data set, supplied as [pk()] argument `data`
#' - `data_settings`: Instructions for data pre-processing (a named list of arguments to [preprocess_data()]), supplied in [pk()] arguments `mapping` and `data_settings`
#' - `scales`: Instructions for data scaling and/or transformation. A list with elements named `conc` and `time`, where each element contains the scaling/transformation to apply to the corresponding variable. See [scale_conc()] and [scale_time()].
#' - `optimx_settings`: Instructions for the numerical optimizer: a named list of arguments to [optimx::optimx()]. See [settings_optimx()].
#' - `status`: What stage of the analysis has been applied to this object so far? Options are 1 (meaning the workflow has been set up), 2 (meaning data has been pre-processed), 3 (meaning that pre-fitting is complete), or 4 (meaning that fitting is complete).
#'
#'No data processing, model fitting, or any other analysis is done until you
#'explicitly request it. Until then, the `pk` object remains just a set of data and
#'instructions. This allows you to specify the instructions for each analysis
#'step without regard for the actual order of the analysis steps, and to
#'overwrite previous instructions, without having to re-do the model fitting
#'each time you add or change a set of instructions. This is particularly useful
#'if you are working in interactive mode at the R console.
#'
#'For example, you might write at the console
#'
#'`my_pk <- pk(my_data) + stat_model(model = "1comp") + data_settings(impute_loq
#'= TRUE)`
#'
#'This is OK even though `data_settings` provides instructions for data
#'pre-processing, a step that comes *before* model fitting. Internally, the `pk`
#'object will put the instructions in the right order.
#'
#'You might then realize that you also want to fit a 2-compartment model to the
#'same data set. You can simply write
#'
#'`my_pk <- my_pk + stat_model(model = "2comp")`
#'
#'Then you might realize that you actually wanted to dose-normalize the
#'concentration data before fitting the models. You can do that simply by
#'writing
#'
#'`my_pk <- my_pk + scale_conc(normalize = "dose")`
#'
#'
#'Now, you are pretty sure that is the final set of instructions. You can
#'actually do the fit as follows:
#'
#'`my_pk <- fit(my_pk)`
#'
#'Now, the following steps will occur:
#'
#' - Data pre-processing, as instructed by [data_settings()]
#' - Data scaling and transformation,as instructed by [scale_conc()] and/or [scale_time()]: Here, concentrations will be dose-normalized
#' - Model pre-fitting
#'      - Automatic determination of whether to fit oral model, IV model, or both, depending on whether oral and IV data are available.
#'      - Automatic checks on whether data are sufficient to proceed with model fitting (e.g., are there more observations than parameters to be estimated?)
#'      - Automatic determination of parameter bounds
#'      - Automatic determination of parameter starting guesses
#' - Model fitting, as instructed by [stat_model()]
#'     - Numerical optimization by maximizing the log-likelihood function
#'     - Calculation of uncertainty in the optimized parameter values using an approximation to the Hessian (the matrix of second derivatives)
#'     - Calculation of Akaike Information Criterion and Bayesian Information Criterion for the fitted model
#' - Model comparison (if more than one model has been fitted), by choosing the model with the lowest AIC value
#'
#'`my_pk` will be modified to contain the results of each of these steps:
#'
#' - the pre-processed, scaled/transformed data (in `my_pk$data`)
#' - for each fitted model: (named elements in `my_pk$models`)
#'     - the parameter bounds and starting guesses (in `my_pk$model$[model name]$parDF`)
#'     - a named vector of the estimated coefficients (in `my_pk$model$[model name]$coefficients`)
#'     - A numerical vector of the residuals for the fitted model (observed - predicted concentrations). If any scaling/transformation was applied, the residuals will be in the transformed units.
#'     - AIC and BIC (in `my_pk$model$[model name]$AIC` and `my_pk$model$[model name]$BIC`))
#' - the winning model by lowest AIC (if more than one model was fitted)
#'
#'You may do these steps one at a time if you wish, using the following methods:
#'
#' - Data pre-processing, including scaling/transformation: [preprocess.pk()]
#' - Model pre-fitting: [prefit.pk()]
#' - Model fitting: [fit.pk()]
#' - Model comparison: [model_compare.pk()]
#'
#'The `pk` object
#'
#'# Mappings
#'
#'Your input data can have any variable names you like. However, internally,
#'`invivopkfit` needs to use a set of "standard", harmonized variable names
#'(e.g., it refers to the variable containing measured tissue concentrations as
#'`Conc`; the variable containing observed time points as `Time`; and the
#'variable containing administered doses as `Dose`). In effect, `invivopkfit`
#'needs to rename the input data, and produce a new `data.frame` that uses these
#'internal harmonized variable names.
#'
#'In order to know which variable names in the input data correspond to each of
#'the internal harmonized variable names, we need to set up a mapping between
#'the internal harmonized variable names and the original variable names.
#'
#'The simplest, most flexible way to set up this mapping is by (ab)using a call
#'to [ggplot2::aes()]. In the context of [ggplot2::ggplot2-package()], you would
#'use [ggplot2::aes()] to set up mappings to `ggplot2`'s "aesthetics", internal
#'harmonized variable names which it uses for plotting: *e.g.*, `x`, `y`,
#'`color`, `size`, `shape`, and so on. In the context of
#'[invivopkfit-package()], we are setting up mappings to `invivopkfit`'s
#'internal harmonized variable names which it uses in model fitting. These
#'"`invivopkfit` aesthetic" variables are as follows:
#'
#'
#' - `Chemical`: A `character` variable containing the chemical identifier. All rows of `data` should have the same value for `Chemical`.
#' - `Species`: A `character` variable containing the name of the species for which the data was measured.  All rows of `data` should have the same value for `Species`.
#' - `Reference`: A `character` variable containing a unique identifier for the data reference (e.g., a single publication).
#' - `Subject`: A `character` variable containing a unique identifier for the subject associated with each observation (an individual animal or group of animals).
#' - `N_Subjects`: A `numeric` variable; an integer giving the number of individual animals represented by this observation. (Some data sets report tissue concentrations for individual animals, in which case `N_Subjects` will be 1; others report average tissue concentrations for groups of multiple animals, in which case `N_Subjects` will be greater than 1.)
#' - `Weight`: A `numeric` variable giving the subject's body weight.
#' - `Weight.Units`: A `character` variable giving the units of body weight.
#' - `Route`: A `character` variable denoting the route of administration. Either `po` (oral) or `iv` (intravenous). Other routes are not currently supported.
#' - `Dose`: A `numeric` variable giving the dose administered.
#' - `Dose.Units`: A `character` variable giving the units of the administered doses.
#' - `Time`: A `numeric` variable giving the time of tissue collection.
#' - `Time.Units`: A `numeric` variable giving the units of `Time`.
#' - `Media`: A `character` variable giving the tissue that was analyzed. Either `blood` or `plasma`. Other tissues are not currently supported.
#' - `Value`: A `numeric` variable giving the tissue concentration in units of mg/L. If `N_Subjects > 1`, `Value` is assumed to represent the mean tissue concentration for this group of subjects. If the tissue concentration was below the limit of quantification (LOQ), this value may be `NA_real_`.
#' - `Value_SD`: A `numeric` variable giving the standard deviation of the tissue concentration in units of mg/L, if available and relevant. If `N_Subjects > 1`, `Value_SD` is assumed to represent the standard deviation of tissue concentrations for this group of subjects. If `N_Subjects == 1`, then `Value_SD` may be `NA_real_`.
#' - `LOQ`: A `numeric` variable giving the limit of quantification applicable to this tissue concentration in units of mg/L, if available.
#' - `Value.Units`: A `character` variable giving the units of `Value`, `Value_SD`, and `LOQ`.
#'
#'You may additionally include mappings to other variable names of your choice, which will appear in
#'the `pk` object in `pk$data` after the analysis is done. The following
#'variable names are reserved for internal use (i.e., they are automatically assigned by [preprocess_data()]:
#'
#' - `Conc`: This is assigned as the greater of `Value` and `LOQ`.
#' - `Detect`: This is a logical variable, `TRUE` if `Conc > LOQ` and `FALSE` otherwise.
#' - `iv`: This is a logical variable, `TRUE` if `Route %in% 'iv'` and `FALSE` otherwise.
#'
#'As with usual calls to [ggplot2::aes()], you should provide the variable names
#'without quoting them. For example, use `ggplot2::aes(Chemical = my_chem)`. Do
#'*not* use `ggplot2::aes("Chemical" = "my_chem")`.
#'
#'Also, as with usual calls to [ggplot2::aes()], you may also specify that any
#'of the "`invivopkfit` aesthetic" variables should be mapped to a constant
#'value, rather than to a variable in `data`. For example, imagine that you
#'don't have a column in `data` that encodes the units of body weight, but you
#'know that all body weights are provided in units of kilograms. You could
#'specify `mapping = ggplot2::aes(Chemical = my_dtxsid, Species = my_species,
#'Weight = my_weight, Weight.Units = "kg")` to map `Weight.Units` to a fixed
#'value of "kg".
#'
#'Finally, as with usual calls to [ggplot2::aes()], you may specify mappings as
#'expressions that use variable names in `data`. For example, if the
#'species-name variable in `data` sometimes says "rat", sometimes "Rat",
#'sometimes "RAT", you might want to harmonize the capitalization. You can do
#'that easily by specifying `mapping = ggplot2::aes(Chemical = my_dtxsid, Species
#'= tolower(my_species)`.
#'
#'# Data
#'
#'`data` should contain data for only one `Chemical` and one `Species`. It may
#'contain data for multiple `Route` and/or `Media` values, which can be fitted
#'simultaneously. However, `Route` values should be either `"po"` (oral bolus
#'administration) or `"iv"` (IV bolus administration), and `Media` values should be
#'either `"blood"` or `"plasma"`.
#'
#'# Data settings
#'
#'The default value of argument `data_settings` is a named list with the
#'following elements:
#'- `ratio_conc_to_dose` Ratio between the mass units used to report the
#'concentration data and the mass units used to report the dose. Default 1. For
#'example, concentration reported in ug/L and dose reported in mg/kg/day would
#'require `ratio_conc_to_dose = 0.001`, because 1 ug/1 mg = 1e-6 g / 1e-3 g =
#'0.001.
#' - `routes_keep` List of routes to retain. Default `c("po", "iv")` to
#'retain only oral and IV administration data.
#'- `media_keep` List of media to retain. Default `c("blood", "plasma")` to
#'retain only concentrations in blood and plasma.
#' - `impute_loq` Logical: TRUE to impute values for missing LOQs; FALSE to
#'leave them alone.
#' - `impute_sd` Logical: TRUE to impute values for missing sample SDs for
#'multi-subject observations; FALSE to leave them alone
#' - `suppress.messages` Logical: Whether to suppress verbose messages.
#'Default FALSE, to be verbose.
#'
#'# Scales
#'
#'The optional argument `scales` is a named list with two elements:
#'
#' - `conc` is itself a list with two elements, `normalize` and `trans`. By default, both are `"identity"`. See [scale_conc()] for options.
#' - `time` is itself a list with one element, `trans`. By default, it is `"identity"`. See [scale_time()] for options.
#'
#'You can always modify `scales` later by using [scale_conc()] and/or
#'[scale_trans()], and it may actually be easier to do it that way. For example,
#'`pk(data = my_df, scales = list(conc = list(trans = log10)))` does exactly the
#'same thing as `pk(data = my_df) + scale_conc(trans = log10)`.
#'
#'# optimx settings
#'
#'Argument `optimx_settings` is a named list containing settings for the
#'optimizer, with the following elements:
#'
#'  - `method` The method to use, from those implemented in [optimx::optimx()]. Default `"bobyqa"`.
#' - `itnmax` The maximum number of iterations, as in [optimx::optimx()]. Default `1e6`.
#' - `hessian` Whether to compute the Hessian after optimizing, as in [optimx::optimx()]. Default `FALSE`.
#' - `control` A named list of control parameters for the optimizer. See [optimx::optimx()] for options and defaults.  Default here is `list(kkt = FALSE)`.
#'
#'
#'
#'@param data A `data.frame`. The default is an empty data frame.
#'@param mapping A mapping set up by (ab)using [ggplot2::aes()]. Call is of form
#'  `ggplot2::aes(new_variable = old_variable)` `new_variable` represents the
#'  harmonized variable name that will be used within `invivopkfit`;
#'  `old_variable` represents the variable name in `data`. If you want to
#'  provide a fixed/constant value for a `new_variable` rather than taking its
#'  value from a variable in `data`, simply supply that fixed/constant value in
#'  the `old_variable` position.
#'@return An object of class `pk`. The initial `pk` object is a list with
#'  elements `data_orig`, `data_settings`, `scales` and `optimx_settings`.
#'  `data_orig` is the original data set to be fitted, as supplied in the
#'  argument `data`. `data_settings` is a named list containing all the other
#'  input arguments: these provide settings that will be used when the data is
#'  pre-processed before fitting.
#'@author Caroline Ring
#'@export

pk <- function(data = NULL,
               mapping = ggplot2::aes(
                Chemical = Chemical,
                Species = Species,
                Reference = Reference,
                Media = Media,
                Route = Route,
                Dose = Dose,
                Dose.Units = "mg/kg",
                Subject = Subject,
                N_Subjects = N_Subjects,
                Weight = Weight,
                Weight.Units = "kg",
                Time = Time,
                Time.Units = "hours",
                Value = Value,
                Value_SD = value_SD,
                LOQ = LOQ,
                Value.Units = "mg/L"
               )
               ){

  #Check to ensure the mapping contains all required harmonized column names
  mapping_default <- ggplot2::aes(
    Chemical = NA_character_,
    Species = NA_character_,
    Reference = NA_character_,
    Media = NA_character_,
    Route = NA_character_,
    Dose = NA_real_,
    Dose.Units = "mg/kg",
    Subject = NA_character_,
    N_Subjects = NA_real_,
    Weight = NA_real_,
    Weight.Units = "kg",
    Time = NA_real_,
    Time.Units = "hours",
    Value = NA_real_,
    Value_SD = NA_real_,
    LOQ = NA_real_,
    Value.Units = "mg/L"
  )

  missing_aes <- setdiff(names(mapping_default), names(mapping))
  if(!(length(missing_aes)==0)){
    mapping[missing_aes] <- mapping_default[missing_aes]
    warning(paste("'mapping' is missing the following required harmonized variables:",
            paste(missing_aes,
                  collapse = "\n"),
            "These missing required variables will be added according to the following mapping:",
            ggplot2:::print.uneval(mapping_default[missing_aes]),
            sep = "\n"))
  }


#Create the initial pk object
  obj <- list("data_original" = data,
              "mapping" = mapping,
              "status" = 1L
              )

#nd assign it class pk
  class(obj) <-append(class(obj), "pk")

# Add default data settings
  obj <- obj + settings_data()

  #Add default scalings for conc and time
  obj <- obj + scale_conc() + scale_time()

  #Add default stats
  #obj <-obj + stat_model()

  #Add default error model
  obj <- obj + stat_error_model()

  #Add default optimx settings
  obj <- obj + settings_optimx()

#return the initialized pk object
  return(obj)

}

#' Check whether an object is of class `pk`
#'
#' @param obj The object whose class is to be tested
#' @return TRUE if the object conforms to expectations for class
#'   `pk`, FALSE if it does not
#' @export
#' @author Caroline Ring
is.pk <- function(obj){
  return(inherits(obj, "pk"))
}


#' Add a pkproto object to a pk object
#'
#' @param e1 A pk pbject
#' @param e2 A pkproto object
#' @return The pk object, modified by adding the pkproto object
#'@export
#' @author Caroline Ring
"+.pk" <- function(e1, e2) {
  if (missing(e2)) {
    cli::cli_abort(c(
      "Cannot use {.code +} with a single argument",
      "i" = "Did you accidentally put {.code +} on a new line?"
    ))
  }

  # Get the name of what was passed in as e2, and pass along so that it
  # can be displayed in error messages
  e2name <- deparse(substitute(e2))

  if      (is.pk(e1))  add_pk(e1, e2, e2name)
  else if (is.pkproto(e1)) {
    cli::cli_abort(c(
      "Cannot add {.cls pkproto} objects together",
      "i" = "Did you forget to add this object to a {.cls pk} object?"
    ))
  }
}

#' Add various `pkproto` objects to a `pk` object
#'
#'@param pk_obj The `pk` object
#'@param object The `pkproto` object to be added
#'@param objectname The name of the `pkproto` object to be added
#'
#'@return The `pk` object modified by the addition.
#' @export
add_pk <- function(pk_obj, object, objectname) {
  if (is.null(object)) return(pk_obj)

  p <- pk_add(object, pk_obj, objectname)
  p
}

#' Add a `pk_scales` object to a `pk` object.
#'
#' @param object The `pk_scales` object to be added.
#' @param pk_obj The `pk` object to which the `pk_scales` object will be added.
#' @param objectname The name of the `pk_scales` object.
#'
#' @return The `pk` object, modified by the `pk_scales` object.
#' @author Caroline Ring
#' @export
pk_add.pk_scales <- function(object, pk_obj, objectname){
pk_obj$scales[[object$name]] <- object$value
if(pk_obj$status > 1L){
  #with new scaling, everything will change starting from data pre-processing
  pk_obj$status <- 1L
}

return(pk_obj)
}

#' Add a `pk_data_settings` object.
#'
#' @param object The `pk_data_settings` object to be added.
#' @param pk_obj The `pk` object to which the `pk_data_settings` object will be added.
#' @param objectname The name of the `pk_data_settings` object.
#'
#' @return The `pk` object, modified by the `pk_data_settings` object.
#' @author Caroline Ring
#' @export
pk_add.pk_data_settings <- function(object, pk_obj, objectname){

  #New data_settings will *replace* existing ones
  if(!is.null(pk_obj$data_settings)){
    message(paste0(objectname,
                   ": data_settings already present; new data_settings will replace the existing one")
    )
  }

pk_obj$data_settings <- object
if(pk_obj$status > 1L){
  #with new data pre-processing settings, everything will change starting from
  #data pre-processing
  pk_obj$status <- 1L
}
return(pk_obj)
}

#' Add a `pk_optimx_settings` object.
#'
#' @param object The `pk_optimx_settings` object to be added.
#' @param pk_obj The `pk` object to which the `pk_optimx_settings` object will be added.
#' @param objectname The name of the `pk_optimx_settings` object.
#'
#' @return The `pk` object, modified by adding the settings.
#' @author Caroline Ring
#' @export
pk_add.pk_optimx_settings <- function(object, pk_obj, objectname){

  #New optimx_settings will *replace* existing ones
  if(!is.null(pk_obj$optimx_settings)){
    message(paste0(objectname,
                   ": optimx_settings already present; new optimx_settings will replace the existing one")
    )
  }

  pk_obj$optimx_settings <- object
  if(pk_obj$status > 3L){
    #with new optimizer settings, data pre=processing and model pre-fitting
    #should not change, but model fitting will change
    message(paste0(objectname,
                   ": New optimx_settings resets status to level 2 (data preprocessing complete); model pre-fit (level 3) and model fit (level 4) will need to be re-done")
    )
    pk_obj$status <- 3L
  }
  return(pk_obj)
}


#' Add a `pk_stat_model` object.
#'
#' @param object The `pk_stat_model` object to be added.
#' @param pk_obj The `pk` object to which the `pk_stat_model` object will be added.
#' @param objectname The name of the `pk_stat_model` object.
#'
#' @return The `pk` object, modified by adding the `stat_model`.
#' @author Caroline Ring
#' @export
pk_add.pk_stat_model <- function(object, pk_obj, objectname){
  #New stat_models will *replace* existing ones by the same name
  for(this_model in names(object)){
    if(!is.null(pk_obj$stat_model[[this_model]])){
      message(paste0(objectname,
                     ": stat_model for",
                     this_model,
                     "already present; new stat_model will replace the existing one")
      )
    }
    pk_obj$stat_model[[this_model]] <- object[[this_model]]
  }
  if(pk_obj$status > 2L){
    message(paste0(objectname,
                   ": New stat_model resets status to level 2 (data preprocessing complete);\nmodel pre-fit (level 3)  (prefit()) and model fit (level 4) (fit()) will need to be re-done")
    )
    pk_obj$status <- 2L #data pre-processing won't change with addition of new model, but model pre-fit and fit will change
  }

  return(pk_obj)
}

#' Add a `pk_stat_error_model` object.
#'
#' @param object The `pk_stat_error_model` object to be added.
#' @param pk_obj The `pk` object to which the `pk_stat_error_model` object will be added.
#' @param objectname The name of the `pk_stat_error_model` object.
#'
#' @return The `pk` object, modified by adding the `stat_error_model`.
#' @author Caroline Ring
#' @export
pk_add.pk_stat_error_model <- function(object, pk_obj, objectname){

  if(!is.null(pk_obj$stat_error_model)){
    message(paste0(objectname,
                   ": stat_error_model already present; new stat_error_model will replace the existing one")
    )
  }

    pk_obj$stat_error_model <- object

    #data pre-processing won't change with addition of a new error model, but
    #model pre-fit and fit will change
  if(pk_obj$status > 2L){
    message(paste0(objectname,
                    ": New stat_error_model resets status to level 2 (data preprocessing complete); model pre-fit (level 3) and model fit (level 4) will need to be re-done")
    )
    pk_obj$status <- 2L
  }

  return(pk_obj)
}


#' Print a PK object
#'
#' Print a
#'
#' A `pk` object is just a list of data and fitting options. In order to
#' actually perform the optimization and fit the model, you need to call one of
#' the methods to do that -- including [print.pk()], [summarize.pk()],
#' [fit.pk()]. If you just type in a set of instructions like `pk(data =
#' my_data) + stat_model(model = c("flat", "1comp", "2comp")` and hit
#' Enter/Return, then by default R will call the [print.pk()] method. (This is
#' true no matter what you type at the R command line and hit enter -- R will
#' call the appropriate `print` method for an object of that class, or
#' [print.default()] if it can't find a class-specific print method.) Therefore,
#' [print_pk()] does the following:
#'
#' - Pre-processes the data
#' - Does initial data checking and summary (e.g., number of observations by route, media, detect/nondetect)
#' - Determines parameters to be optimized for each specified model, based on the data
#' - Checks data to see whether the selected parameters may be identifiable (e.g., are there more observations than there are parameters to be estimated?)
#' - Determines bounds and starting values for each parameter to be optimized
#' - Performs optimization to estimate parameter values and uncertainties
#' - Adds the optimization results to the `pk` object
#' - Prints the `pk` object
#' - Returns the `pk` object invisibly
#'
#' @param obj A `pk` object
#' @return Invisibly: The `pk` object with added elements containing the optimization results
#' @author Caroline Ring
#' @export
print.pk <- function(obj){
str(obj)
}

#' Pre-process data
#'
#' S3 method for `pk` objects to preprocess data
#'
#' This is the `preprocess_data` S3 method for objects of class `pk`.
#'
#' @param obj A `pk` object
#' @return The same `pk` object, with added elements `data` (containing the
#'   cleaned, gap-filled data) and `data_info` (containing summary information
#'   about the data, e.g. number of observations by route, media,
#'   detect/nondetect; empirical tmax, time of peak concentration for oral data;
#'   number of observations before and after empirical tmax)
#' @author Caroline Ring
#' @importFrom magrittr `%>%`
#' @export
preprocess_data.pk <- function(obj){

  if(is.null(obj$data_original)){
    message("Original data is NULL")
    obj$data <- NULL
    obj$data_info <- NULL
    obj$status <- 2
    return(obj)
  }else{

    #coerce to data.frame (in case it is a tibble or data.table or whatever)
    data_original <- as.data.frame(obj$data_original)

    #rename variables using the mapping
    data <- as.data.frame(sapply(obj$mapping,
                                 function(x) rlang::eval_tidy(x, data_original),
                                 simplify = FALSE,
                                 USE.NAMES = TRUE)
    )

    #Check to make sure the data include only one Chemical and Species. Stop
    #with an error otherwise.
    chems <- unique(data$Chemical)
    species <- unique(data$Species)

    nchem <- length(chems)
    nspecies <- length(species)

    if(!(nchem %in% 1 & nspecies %in% 1)){
      stop(paste("preprocess_data.pk(): data contains multiple chemicals and/or multiple species.",
           "Unique chemicals in this data:",
           paste(chems, collapse = "; "),
           "Unique Species in this data:",
           paste(species, collapse = "; "),
           sep = "\n"))
    }

    #Check to make sure the data include only route_keep and media_keep
    routes <- unique(data$Route)
    media <- unique(data$Media)

    if( !(all(routes %in% obj$data_settings$routes_keep) &
          all(media %in% obj$data_settings$media_keep)) ){
      stop(paste("preprocess_data.pk(): data contains unsupported media and/or routes.",
                 paste("Supported media:", paste(obj$data_settings$media_keep, collapse = "; ")),
                 paste("Media in data:", paste(media, collapse = "; ")),
                 paste("Supported routes:", paste(obj$data_settings$routes_keep, collapse = "; ")),
                 paste("Routes in data:", paste(routes, collapse = "; ")),
                 sep = "\n"))
    }


    #Check to make sure there is only one value for each set of units
    time_units <- unique(data$Time.Units)
    value_units <- unique(data$Value.Units)
    weight_units <- unique(data$Weight.Units)
    dose_units <- unique(data$Dose.Units)

    if(any(sapply(list(time_units,
                       value_units,
                       weight_units,
                       dose_units),
                  function(x) length(x) > 1)
    )){
      stop(paste("preprocess_data.pk(): data contains multiple units for one or more variables.",
                 paste("Time units:", paste(time_units, collapse = "; "), sep = " ",),
                 paste("Concentration units:", paste(value_units, collapse = "; "), sep = " "),
                 paste("Weight units:", paste(weight_units, collapse = "; "), sep = " "),
                 paste("Dose units:", paste(dose_units, collapse = "; "), sep = " "),
                 sep = "\n"))
    }

    # If data has passed all these initial checks, then proceed with pre-processing

    if(!obj$data_settings$suppress.messages){
      ### display messages describing loaded data
      message(
        paste(
          paste(nrow(data),
                "concentration vs. time observations loaded."),
          paste("Chemicals:", chems),
          paste("Species:", species),
          paste("Routes:", routes),
          paste("Media:", media),
          sep = "\n"
        )
      )
    }

    ### Coerce all 'Value' values to be numeric and say so
    if (!is.numeric(data$Value))
    {
      value_num <- as.numeric(data$Value)
      old_na <- sum(is.na(data$Value) | !nzchar(data$Value))
      new_na <- sum(is.na(value_num))
      if(!obj$data_settings$suppress.messages){
        message(paste0("Column \"Value\" converted from ",
                       class(data$Value),
                       " to numeric. ",
                       "Pre-conversion NAs and blanks: ",
                       old_na,
                       ". Post-conversion NAs: ",
                       new_na, "."))
      }
      data$Value <- value_num
      rm(value_num, old_na, new_na)
    }

    ### coerce 'Dose' values to numeric and say so
    if (!is.numeric(data$Dose))
    {
      dose_num <- as.numeric(data$Dose)
      old_na <- sum(is.na(data$Dose) | !nzchar(data$Dose))
      new_na <- sum(is.na(dose_num))
      if(!obj$data_settings$suppress.messages){
        message(paste0("Column \"Dose\" converted from ",
                       class(data$Dose),
                       " to numeric. ",
                       "Pre-conversion NAs and blanks: ",
                       old_na,
                       ". Post-conversion NAs: ",
                       new_na, "."))
      }
      data$Dose <- dose_num
      rm(dose_num, old_na, new_na)
    }

    ### coerce 'Time' values to numeric and say so
    if (!is.numeric(data$Time))
    {
      time_num <- as.numeric(data$Time)
      old_na <- sum(is.na(data$Time) | !nzchar(data$Time))
      new_na <- sum(is.na(time_num))
      if(!obj$data_settings$suppress.messages){
        message(paste0("Column \"Time\" converted from ",
                       class(data$TIme),
                       " to numeric. ",
                       "Pre-conversion NAs and blanks: ",
                       old_na,
                       ". Post-conversion NAs: ",
                       new_na, "."))
      }
      data$Time <- time_num
      rm(time_num, old_na, new_na)
    }

    ### Coerce Species, Route, and Media to lowercase
    data$Species <- tolower(data$Species)
    data$Route <- tolower(data$Route)
    data$Media <- tolower(data$Media)

    #Coerce any non-positive Value, LOQ, or Value_SD to NA
    data[(data$Value<=0) %in% TRUE, "Value"] <- NA_real_
    #Coerce any non-positive LOQ to NA
    data[(data$LOQ<=0) %in% TRUE, "LOQ"] <- NA_real_
    #Coerce any non-positive Value_SD to NA
    data[(data$Value_SD<=0) %in% TRUE, "Value_SD"] <- NA_real_

    # Impute LOQ
    data$LOQ_orig <- data$LOQ
    if(obj$data_settings$impute_loq %in% TRUE){
      if(any(is.na(data$LOQ))){
        if(!obj$data_settings$suppress.messages){
          message(paste0("Estimating missing LOQs as ",
                         obj$data_settings$calc_loq_factor,
                         "* minimum detected Value for each unique combination of ",
                         paste(
                           sapply(obj$data_settings$loq_group,
                               as_label),
                           collapse = " + ")
                         )
          )
        }
        data <- do.call(dplyr::group_by,
                        args =c(list(data),
                                obj$data_settings$loq_group)) %>%
          dplyr::mutate(LOQ_orig = LOQ,
                 LOQ = dplyr::if_else(is.na(LOQ_orig),
                               min(LOQ_orig, na.rm = TRUE) *
                                 obj$data_settings$calc_loq_factor,
                               LOQ_orig)
                 ) %>%
          dplyr::ungroup() %>%
          as.data.frame()
      }

    } #end if impute_loq %in% TRUE

    if(!obj$data_settings$suppress.messages){
      message(paste0("Converting 'Value' values of less than LOQ to NA.\n",
                     sum((data$Value <  data$LOQ) %in% TRUE),
                     " values will be converted."))
    }
    data[(data$Value < data$LOQ) %in% TRUE, "Value"] <- NA_real_

    #Remove any remaining cases where both Value and LOQ are NA
    if(any(is.na(data$Value) & is.na(data$LOQ))){
      if(!obj$data_settings$suppress.messages){
        message(paste0("Removing observations where both Value and LOQ were NA.\n",
                       sum(is.na(data$Value) & is.na(data$LOQ)),
                       " observations will be removed."))
      }
      data <- subset(data,
                     !(is.na(data$Value) & is.na(data$LOQ)))

      if(!obj$data_settings$suppress.messages){
        message(paste(dim(data)[1], "observations of",
                      length(unique(data$Chemical)), "unique chemicals,",
                      length(unique(data$Species)), "unique species, and",
                      length(unique(data$Reference)), "unique references remain."))
      }
    }

    # For any cases where N_Subjects is NA, impute N_Subjects = 1
    if(any(is.na(data$N_Subjects))){
      data$N_Subjects_orig <- data$N_Subjects
      if(!obj$data_settings$suppress.messages){
        message(
          paste0(
            "N_Subjects is NA for ",
            sum(is.na(data$N_Subjects)),
            " observations. It will be assumed = 1."
          )
        )
      }

      data[is.na(data$N_Subjects), "N_Subjects"] <- 1
    }

    #for anything with N_Subjects == 1, set Value_SD to 0
    data[data$N_Subjects == 1, "Value_SD"] <- 0

    # Impute missing SDs
    data$Value_SD_orig <- data$Value_SD
    if(obj$data_settings$impute_sd %in% TRUE){
      if(any((data$N_Subjects >1) %in% TRUE & is.na(data$Value_SD))){
        if(!obj$data_settings$suppress.messages){
          #number of SDs to be estimated
          n_sd_est <- sum(
            (data$N_Subjects >1) %in% TRUE &
              is.na(data$Value_SD)
          )
          message(paste0("Estimating missing concentration SDs for multi-subject data points as ",
                         "minimum non-missing SD for each unique combination of variables in",
                         paste(
                           sapply(obj$data_settings$sd_group,
                                  as_label),
                           collapse = " + "),
                         ". If all SDs are missing for such a unique combination, ",
                         "SD will be imputed equal to mean. ",
                         n_sd_est, " missing SDs will be estimated."))
        }
        data <- do.call(dplyr::group_by,
                        args = c(list(data),
                                 obj$data_settings$sd_group)) %>%
          dplyr::mutate(Value_SD_orig = Value_SD,
                 Value_SD = dplyr::if_else(is.na(Value_SD_orig) &
                                      N_Subjects > 1,
                                    dplyr::if_else(all(is.na(Value_SD_orig)),
                                            mean(Value, na.rm = TRUE),
                                            min(Value_SD_orig, na.rm = TRUE)),
                               Value_SD_orig)
          ) %>% dplyr::ungroup() %>%
          as.data.frame()
      }

    }#end if impute_sd %in% TRUE

    #Remove any remaining multi-subject observations where SD is NA
    #(with imputing SD = Mean as a fallback, this will only be cases where Value was NA)
    if(any((data$N_Subjects >1) %in% TRUE & is.na(data$Value_SD))){
      if(!obj$data_settings$suppress.messages){
        message(paste0("Removing observations with N_Subjects > 1 where reported SD is NA.\n",
                       sum((data$N_Subjects >1) %in% TRUE &
                             is.na(data$Value_SD)),
                       " observations will be removed."))
      }
      data <- subset(data,
                     !((N_Subjects >1) %in% TRUE & is.na(Value_SD))
      )

      if(!obj$data_settings$suppress.messages){
        message(paste(dim(data)[1], "observations of",
                      length(unique(data$Chemical)), "unique chemicals,",
                      length(unique(data$Species)), "unique species, and",
                      length(unique(data$Reference)), "unique references remain."))
      }
    }

    #Remove any remaining multi-subject observations where Value is NA
    if(any((data$N_Subjects >1) %in% TRUE & is.na(data$Value))){
      if(!obj$data_settings$suppress.messages){
        message(paste0("Removing observations with N_Subjects > 1 where reported Value is NA.\n",
                       sum((data$N_Subjects >1) %in% TRUE &
                             is.na(data$Value)),
                       " observations will be removed."))
      }
      data <- subset(data,
                     !((N_Subjects >1) %in% TRUE & is.na(Value))
      )

      if(!obj$data_settings$suppress.messages){
        message(paste(dim(data)[1], "observations of",
                      length(unique(data$Chemical)), "unique chemicals,",
                      length(unique(data$Species)), "unique species, and",
                      length(unique(data$Reference)), "unique references remain."))
      }
    }

    #Remove any NA time values
    if(any(is.na(data$Time))){
      if(!obj$data_settings$suppress.messages){
        message(paste0("Removing observations with NA time values.\n",
                       sum(is.na(data$Time)),
                       " observations will be removed."))
      }
      data <- subset(data, !is.na(Time))

      if(!obj$data_settings$suppress.messages){
        message(paste(dim(data)[1], "observations of",
                      length(unique(data$Reference)), "unique references remain."))
      }
    }

    #Remove any Dose = 0 observations
    if(any(data$Dose <= .Machine$double.eps)){
      if(!obj$data_settings$suppress.messages){
        message(paste0("Removing observations with Dose == 0.\n",
                       sum(data$Dose <= .Machine$double.eps),
                       " observations will be removed."))
      }
      data <- subset(data, data$Dose > .Machine$double.eps)

      if(!obj$data_settings$suppress.messages){
        message(paste(dim(data)[1], "observations of",
                      length(unique(data$Reference)), "unique references remain."))
      }
    }

  #apply time transformation
  data$Time_orig <- data$Time
  data$Time.Units_orig <- data$Time.Units

  #first, default to identity transformation if none is specified
  if(is.null(obj$scales$time$new_units)){
    obj$scales$time$new_units <- "identity"
  }


  from_units <- unique(data$Time.Units)
  to_units <- ifelse(obj$scales$time$new_units %in% "identity",
                     from_units,
                     obj$scales$time$new_units)



  if(obj$scales$time$new_units %in% "auto"){
    to_units <- auto_units(y = data$Time,
                     from = from_units)
  }

  if(!obj$data_settings$suppress.messages){
    message(paste("Converting time from",
                  from_units,
                  "to",
                  to_units))
  }
  data$Time_trans <- tryCatch(convert_time(x = data$Time,
                                     from = from_units,
                                     to = to_units,
                                     inverse = FALSE),
                        error = function(err){
                          warning(paste("invivopkfit::preprocess_data.pk():",
                                  "Error in transforming time using convert_time():",
                                  err$message))
                          return(NA_real_)
                        })
  data$Time_trans.Units <- to_units

  #If Conc and Detect do not already exist, create them
  if(!("Conc" %in% names(data))){
    if(!obj$data_settings$suppress.messages){
      message(paste("Harmonized variable `Conc` does not already exist.",
                    "It will be created as",
                    "`pmax(Value, LOQ, na.rm = TRUE)`."))
    }
    data$Conc <- pmax(data$Value, data$LOQ, na.rm = TRUE)
  }

  if(!("Detect" %in% names(data))){
    if(!obj$data_settings$suppress.messages){
      message(paste("Harmonized variable `Detect` does not already exist.",
                    "It will be created as",
                    "`is.na(Value)`."))
    }
    data$Detect <- !is.na(data$Value)
  }

  if(!("Conc_SD" %in% names(data))){
    if(!obj$data_settings$suppress.messages){
      message(paste("Harmonized variable `Conc_SD` does not already exist.",
                    "It will be created as",
                    "`ifelse(data$Detect,
                        data$Value_SD,
                        NA_real_)"))
    }
    data$Conc_SD <- data$Value_SD
  }

  data$Conc.Units <- data$Value.Units

  #apply concentration transformation.
  #
  #.conc is a placeholder that refers to any concentration variable (Conc,
  #value, Value_SD, Conc_SD, LOQ).

  #apply conc transformation
  #use tidy evaluation: Specify an expression to evaluate in the context of a data.frame
  #Transform Conc (this is either the measured value or the LOQ, depending on Detect)
  #If no conc transformation specified, assume identity
  if(is.null(obj$scales$conc$expr)){
    obj$scales$conc$ratio_conc_dose <- 1
    obj$scales$conc$dose_norm <- FALSE
    obj$scales$conc$log10_trans <- FALSE
    obj$scales$conc$expr <- rlang::new_quosure(quote(.conc),
                                               env = caller_env())
  }

  if(!obj$data_settings$suppress.messages){
    message(paste("Applying transformations to concentration variables:",
                  rlang::as_label(obj$scales$conc$expr)))
  }

  data$Conc_trans <- rlang::eval_tidy(expr = obj$scales$conc$expr,
                         data = cbind(data,
                                      data.frame(.conc = data$Conc)
                                      ))
  #Transform Conc_SD
  data$Conc_SD_trans <- rlang::eval_tidy(obj$scales$conc$expr,
                                data = cbind(data,
                                             data.frame(.conc = data$Conc_SD)
                                ))
  #Record new conc units
  data$Conc_trans.Units <- gsub(x = rlang::as_label(obj$scales$conc$expr),
                          pattern = ".conc",
                          replacement = unique(data$Value.Units),
                          fixed = TRUE)

  #get the summary data info
  #unique Chemical, Species, References, Studies, Routes
  dat_info <- as.list(unique(data[c("Chemical",
                                    "Species")]))

  dat_info$References_Analyzed <- sort(unique(data$Reference))
  #get a list of studies analyzed
  dat_info$Studies_Analyzed <- sort(unique(data$Study))
  #get a list of routes analyzed
  dat_info$Routes_Analyzed <- sort(unique(data$Route))

  #get a list of media analyzed
  dat_info$Media_Analyzed <- sort(unique(data$Media))

  #get the number of detects and non-detects by route and medium
  dat_info$n_dat <- aggregate(x = list(Detect = data$Detect),
                              by = data[c("Route", "Media")],
                              FUN = function(Detect){
                                c("Detect" = sum(Detect %in% TRUE),
                                  "NonDetect" = sum(Detect %in% FALSE))
                              })
  names(dat_info$n_dat) <- gsub(x = names(dat_info$n_dat),
                                pattern = "Detect.",
                                replacement = "",
                                fixed = TRUE)

  #get time of last detected observation
  if(any(data$Detect %in% TRUE)){
  dat_info$last_detect_time <- max(data[data$Detect %in% TRUE, "Time"])
  }else{
    dat_info$last_detect_time <- 0
  }

  #get time of last observation
  dat_info$last_time <- max(data$Time)

  # add data & data info to object
  obj$data <- data
  obj$data_info <- dat_info

  obj$status <- 2 #preprocessing complete

  return(obj)
}
}

#' Do pre-fit calculations and checks
#'
#' @param obj A `pk` object
#' @return The same `pk` object, but with additional elements added to each item
#'   in `models`, containing the results of pre-fit calculations and checks for
#'   each model.
#' @export
#' @author Caroline Ring
prefit.pk <- function(obj){
  #if preprocessing not already done, do it
  if(obj$status < 1){
    obj <- preprocess_data(obj)
  }

  #get the error model obj$stat_error_model, which defines the number of sigmas that will need to be optimized
#count the number of unique combinations of vars in obj$stat_error_model$error_group
unique_groups <- unique(rlang::eval_tidy(expr = obj$stat_error_model$error_group,
                                         data = obj$data))
n_sigma <- nrow(unique_groups)

obj$stat_error_model$n_sigma <- n_sigma

#Assign a factor variable denoting sigma group to each observation in obj$data
#This tells us which sigma applies to which observation
obj$stat_error_model$data_sigma_group <- interaction(
  lapply(
    obj$stat_error_model$error_group,
    function(x){
      rlang::eval_tidy(x, data = obj$data)
    }
  )
)

#get bounds and starting points for each error sigma to be fitted
sigma_DF <- data.frame(param_name = paste("sigma",
                                          levels(obj$stat_error_model$data_sigma_group),
                                          sep = "_"),
                       param_units = unique(obj$data$Conc_trans.Units),
                       optimize_param = TRUE,
                       use_param = TRUE,
                       lower_bound = .Machine$double.eps)

#get upper bound: standard deviation of the transformed concentration
#(Conc_trans) in each group
sigma_DF$upper_bound <- tapply(X = obj$data$Conc_trans,
               INDEX = obj$stat_error_model$data_sigma_group,
               FUN = sd,
               na.rm = TRUE,
               simplify = TRUE)

#get starting value for sigma: say, 0.5 of the upper bound
sigma_DF$start <- 0.5*sigma_DF$upper_bound

#assign rownames to sigma_DF
rownames(sigma_DF) <- sigma_DF$param_name

#assign sigma_DF to the `pk` object
obj$stat_error_model$sigma_DF <- sigma_DF

n_sigma <- nrow(sigma_DF)
  #for each model to be fitted:
  for (this_model in names(obj$stat_model)){
    #get parameters to be optimized, bounds, and starting points
    #by evaluating params_fun for this stat_model
    obj$stat_model[[this_model]]$par_DF <- do.call(obj$stat_model[[this_model]]$params_fun,
                                                   args = c(list(obj$data),
                                                            obj$stat_model[[this_model]]$params_fun_args))
    #check whether there are enough observations to optimize the requested parameters plus sigmas
    #number of parameters to optimize
    n_par <- sum(obj$stat_model[[this_model]]$par_DF$optimize_param)
    #number of detected observations
    n_detect <- sum(obj$data$Detect)
    obj$stat_model[[this_model]]$status <- ifelse(
      n_detect <= (n_par + n_sigma),
      "abort",
     "continue"
    )

    obj$stat_model[[this_model]]$status_reason <- ifelse(
      n_detect <= (n_par + n_sigma),
      paste0("Number of detects (",
             n_detect,
             ") is less than or equal to number of parameters to optimize (",
             n_par, ") plus number of error SDs to optimize (",
             n_sigma,
             ")"),
      paste0("Number of detects (",
             n_detect,
             ") is greater than number of parameters to optimize (",
             n_par, ") plus number of error SDs to optimize (",
             n_sigma,
             ")")
    )

  }

  obj$status <- 3 #prefit complete

  return(obj)

}

#' Fit PK model(s) for a `pk` object
#'
#' @param obj A [pk] object.
#' @return The same [pk] object, with element `fit` added to each `$stat_model`
#'   element, reflecting the fitted parameters for the corresponding model.
#' @export
#' @author Caroline Ring
fit.pk <- function(obj){
suppress.messages <- obj$data_settings$suppress.messages
  #For each model:
  for (this_model in names(obj$stat_model)){

    if(obj$stat_model[[this_model]]$status %in% "continue"){

  #Pull par_DF to get which params to optimize, bounds, and starting values
  par_DF <- obj$stat_model[[this_model]]$par_DF
  #Do the same for sigma_DF (the data frame of error SDs)
  sigma_DF <- obj$stat_error_model$sigma_DF

  #Rowbind par_DF and sigma_DF
  par_DF <- rbind(par_DF,
                  sigma_DF)

  #get params to be optimized with their starting points
  opt_params <- par_DF[par_DF$optimize_param %in% TRUE,
                       "start"]
  names(opt_params) <- par_DF[par_DF$optimize_param %in% TRUE,
                              "param_name"]

  #param lower bounds (only on params to be optimized)
  lower_params <- par_DF[par_DF$optimize_param %in% TRUE,
                         "lower_bound"]
  names(lower_params) <- par_DF[par_DF$optimize_param %in% TRUE,
                                "param_name"]

  #param upper bounds (only on params to be optimized)
  upper_params <- par_DF[par_DF$optimize_param %in% TRUE,
                         "upper_bound"]
  names(upper_params) <- par_DF[par_DF$optimize_param %in% TRUE,
                                "param_name"]

  #get params to be held constant, if any
  if(any(par_DF$optimize_param %in% FALSE &
         par_DF$use_param %in% TRUE)){
    const_params <- par_DF[par_DF$optimize_param %in% FALSE &
                             par_DF$use_param %in% TRUE,
                           "start"]

    names(const_params) <- par_DF[par_DF$optimize_param %in% FALSE &
                                    par_DF$use_param %in% TRUE,
                                  "param_name"]
  }else{
    const_params <- NULL
  }

  fitdata <- obj$data
  #Now call optimx::optimx() and do the fit
  optimx_out <- tryCatch({

    do.call(
      optimx::optimx,
      args = c(
        #
        list(par = opt_params,
             fn = log_likelihood,
             lower = lower_params,
             upper = upper_params),
        #method and control
        obj$optimx_settings,
        #... additional args to log_likelihood
        list(
          const_params = const_params,
          fitdata = fitdata,
          data_sigma_group = obj$stat_error_model$data_sigma_group,
          modelfun = obj$stat_model[[this_model]]$conc_fun,
          scales_conc = obj$scales$conc,
          negative = TRUE,
          force_finite = TRUE
        ) #end list()
      ) #end args = c()
    ) #end do.call

  },
  error = function(err){
    return(paste0("Error from optimx::optimx(): ",
                  err$message))
  })

  #Save the fitting results for this model
  obj$stat_model[[this_model]]$fit <- optimx_out

  #loop over rows of optimx_out (i.e. over optimx methods), if any
  if(!is.null(nrow(optimx_out))){
  #For each row in optimx output, evaluate Hessian and attempt to invert
   for (i in 1:nrow(optimx_out)){
     this_method <- rownames(optimx_out)[i]
     #Get SDs from Hessian evaluated at the fit parameters
     npar <- attr(optimx_out, "npar")

     fit_par <- unlist(optimx_out[i, 1:npar])
     #Calculate Hessian using function from numDeriv
     numhess <- numDeriv::hessian(func = function(x){
       log_likelihood(x,
                      const_params = const_params,
                      fitdata = fitdata,
                      data_sigma_group = obj$stat_error_model$data_sigma_group,
                      modelfun = obj$stat_model[[this_model]]$conc_fun,
                      scales_conc = obj$scales$conc,
                      negative = TRUE,
                      force_finite = TRUE)
     },
     x = fit_par,
     method = 'Richardson')

     #try inverting Hessian to get SDs
     sds <- tryCatch(diag(solve(numhess)) ^ (1/2),
                     error = function(err){
                       #if hessian can't be inverted
                       if (!suppress.messages) {
                         message(paste0("Hessian can't be inverted, ",
                                        "using pseudovariance matrix ",
                                        "to estimate parameter uncertainty."))
                       }
                       #pseudovariance matrix
                       #see http://gking.harvard.edu/files/help.pdf
                       suppressWarnings(tmp <- tryCatch(
                         diag(chol(MASS::ginv(numhess),
                                   pivot = TRUE)) ^ (1/2),
                         error = function(err){
                           if (!suppress.messages) {
                             message(paste0("Pseudovariance matrix failed,",
                                            " returning NAs"))
                           }
                           rep(NA_real_, nrow(numhess))
                         }
                       )
                       )
                       return(tmp)
                     })
     names(sds) <- names(optimx_out)[1:npar]

     obj$stat_model[[this_model]]$fit_hessian[[this_method]] <- numhess
     obj$stat_model[[this_model]]$fit_sds[[this_method]] <- sds

    #  #calculate also: log-likelihood, AIC, BIC
    # obj$stat_model[[this_model]]$fit_AIC[[this_method]] <-
    #   obj$stat_model[[this_model]]$fit_BIC[[this_method]] <-
   } #end loop over rows of optimx_out
  } #end check that there *were* any rows of optimx_out
}else{ #if status for this model was "abort", then abort fit and return NULL
    obj$stat_model[[this_model]]$fit <- "Fit aborted because status %in% 'abort'"
    }
  } #end loop over models

  obj$status <- 4 #fitting complete
  return(obj)
}

#' Get coefficients
#'
#' Extract coefficients from a fitted `pk` object
#'
#' @param obj A [pk] object
#' @param model Optional: Specify one or more of the fitted models whose
#'   coefficients to return. If NULL (the default), coefficients will be returned for all of
#'   the models in `obj$stat_model`.
#' @param method Optional: Specify one or more of the [optimx::optimx()] methods
#'   whose coefficients to return. If NULL (the default), coefficients will be returned for
#'   all of the models in `obj$optimx_settings$method`.
#' @return A named list of numeric matrixes. There is one list element named for
#'   each model in `model`. Each list element is a matrix with as many
#'   rows as items in `method`. The row names are the method names. The matrix column names are
#'   the names of the fitted parameters, including any error standard deviation
#'   hyperparameters (whose names begin with "sigma").
#' @export
#' @author Caroline Ring
coef.pk <- function(obj,
                    model = NULL,
                    method = NULL){
  if(is.null(model)) model <- names(obj$stat_model)
  if(is.null(method)) method <- obj$optimx_settings$method

  sapply(obj$stat_model[model],
         function(this_model){
           npar <- attr(this_model$fit, "npar")
           fit_par <- this_model$fit[method, 1:npar]

           #Add any "constant" params
           if(any(this_model$par_DF$optimize_param %in% FALSE &
                  this_model$par_DF$use_param %in% TRUE)){
           const_parDF <- subset(this_model$par_DF,
                               optimize_param %in% FALSE &
                                 use_param %in% TRUE)[c("param_name",
                                                        "start")]
           const_par <- const_parDF[["start"]]
           names(const_par) <- const_parDF[["param_name"]]
           const_par <- do.call(rbind,
                                replicate(nrow(fit_par),
                                          const_par,
                                          simplify = FALSE))
           fit_par <- as.matrix(cbind(fit_par, const_par))

           #to be implemented: return fit SDs as well as fit values??
           # fit_sd <- do.call(rbind,
           #                   this_model$fit_sds[method])
          #  const_sd <- rep(NA_real_, length(const_par))
          #  names(const_sd) <- names(const_par)
          # fit_sd <- cbind(fit_sd, const_sd)
          #  rownames(fit_sd) <- paste(rownames(fit_sd),
          #                            "fit_sd",
          #                            sep = ".")

           # fit_mat <- rbind(fit_par,
           #                  fit_sd)
           fit_mat <- fit_par
           }
           return(fit_mat)
         },
         USE.NAMES = TRUE,
         simplify = FALSE)
}

#' Get coefficient standard deviations
#'
#' Extract coefficient standard deviations from a fitted `pk` object
#'
#' @param obj A [pk] object
#' @param model Optional: Specify one or more of the fitted models whose
#'   coefficients to return. If NULL (the default), coefficients will be returned for all of
#'   the models in `obj$stat_model`.
#' @param method Optional: Specify one or more of the [optimx::optimx()] methods
#'   whose coefficients to return. If NULL (the default), coefficients will be returned for
#'   all of the models in `obj$optimx_settings$method`.
#' @return A named list of numeric matrixes. There is one list element named for
#'   each model in `model`. Each list element is a matrix with as many
#'   rows as items in `method`. The row names are the method names. The matrix column names are
#'   the names of the fitted parameters, including any error standard deviation
#'   hyperparameters (whose names begin with "sigma").
#' @export
#' @author Caroline Ring
coef_sd.pk <- function(obj,
                    model = NULL,
                    method = NULL){
  if(is.null(model)) model <- names(obj$stat_model)
  if(is.null(method)) method <- obj$optimx_settings$method

  sapply(obj$stat_model[model],
         function(this_model){
           npar <- attr(this_model$fit, "npar")
           fit_par <- this_model$fit[method, 1:npar]

           #Add any "constant" params
           if(any(this_model$par_DF$optimize_param %in% FALSE &
                  this_model$par_DF$use_param %in% TRUE)){
             const_parDF <- subset(this_model$par_DF,
                                   optimize_param %in% FALSE &
                                     use_param %in% TRUE)[c("param_name",
                                                            "start")]
             const_par <- const_parDF[["start"]]
             names(const_par) <- const_parDF[["param_name"]]
             const_par <- do.call(rbind,
                                  replicate(nrow(fit_par),
                                            const_par,
                                            simplify = FALSE))
             fit_par <- as.matrix(cbind(fit_par, const_par))
             #print(const_par)

             fit_sd <- do.call(rbind,
                               this_model$fit_sds[method])

              const_sd <- rep(NA_real_, length(const_par))
              names(const_sd) <- const_parDF[["param_name"]]
              # print(fit_sd)
              # print(const_sd)
              const_sd <- matrix(data = NA_real_,
                                 nrow = nrow(fit_sd),
                                 ncol = ncol(const_par))
              colnames(const_sd) <- colnames(const_par)
             fit_sd <- cbind(fit_sd, const_sd)


             fit_mat <- fit_sd
           }
           return(fit_mat)
         },
         USE.NAMES = TRUE,
         simplify = FALSE)
}

#' Get predictions
#'
#' Extract predictions from a fitted `pk` object.
#'
#' @param obj A [pk] object.
#' @param newdata Optional: A `data.frame` with new data for which to make
#'   predictions. If NULL (the default), then predictions will be made for the
#'   data in `obj$data`. `newdata` is required to contain at least the following
#'   variables: `Time`, `Dose`, `Route`, and `Media`. `Time` will be transformed
#'   according to the transformation in `obj$scales$time` before predictions are
#'   made.
#' @param model Optional: Specify one or more of the fitted models for which to
#'   make predictions. If NULL (the default), predictions will be returned for
#'   all of the models in `obj$stat_model`.
#' @param method Optional: Specify one or more of the [optimx::optimx()] methods
#'   for which to make predictions. If NULL (the default), predictions will be
#'   returned for all of the models in `obj$optimx_settings$method`.
#' @param type Either `"conc"` (the default) or `"auc"`. `type = "conc"`
#'   predicts concentrations; `type = "auc"` predicts area under the
#'   concentration-time curve (AUC).
#' @return A named list of numeric matrixes. There is one list element named for
#'   each model in `obj`'s [stat_model()] element, i.e. each PK model that was
#'   fitted to the data. Each list element is a matrix with the same number of
#'   rows as the data in `obj$data` (corresponding to the rows in `obj$data`),
#'   and as many columns as there were [optimx::optimx()] methods (specified in
#'   [settings_optimx()]). The column names are the method names. Each column
#'   contains the predictions of the model fitted by the corresponding method.
#'   These predictions are concentrations in the same units as
#'   `obj$data$Conc.Units`; any concentration transformations (in
#'   `obj$scale$conc`) are *not* applied.
#' @export
#' @author Caroline Ring
predict.pk <- function(obj,
                       newdata = NULL,
                       model = NULL,
                       method = NULL,
                       type = "conc"
                       ){

  if(is.null(model)) model <- names(obj$stat_model)
  if(is.null(method)) method <- obj$optimx_settings$method

  coefs <- coef(obj = obj,
                model = model,
                method = method)

  if(is.null(newdata)) newdata <- obj$data
  if(!("Time_trans" %in% names(newdata))){
    #transform time if needed
    #first, default to identity transformation if none is specified
    if(is.null(obj$scales$time$new_units)){
      obj$scales$time$new_units <- "identity"
    }

    from_units <- unique(newdata$Time.Units)
    to_units <- ifelse(obj$scales$time$new_units %in% "identity",
                       from_units,
                       obj$scales$time$new_units)


    if(obj$scales$time$new_units %in% "auto"){
      to_units <- auto_units(y = newdata$Time,
                             from = from_units)
    }
    if(!obj$data_settings$suppress.messages){
      message(paste("Converting time from",
                    from_units,
                    "to",
                    to_units))
    }
    newdata$Time_trans <- tryCatch(convert_time(x = newdata$Time,
                                             from = from_units,
                                             to = to_units,
                                             inverse = FALSE),
                                error = function(err){
                                  warning(paste("invivopkfit::predict.pk():",
                                                "Error in transforming time in `newdata` using convert_time():",
                                                err$message))
                                  return(NA_real_)
                                })
  }


  #loop over models
  sapply(model,
         function(this_model){
           this_coef_mat <- coefs[[this_model]]
          apply(this_coef_mat,
                1,
                function(this_coef_row){
                  #get coefficients
                  this_coef <- as.list(this_coef_row)
                  #get model function to be evaluated
                  this_model_fun <- ifelse(type %in% "conc",
                                           obj$stat_model[[this_model]]$conc_fun,
                                           ifelse(type %in% "auc",
                                                  obj$stat_model[[this_model]]$auc_fun,
                                                  NULL)
                  )

                  #evaluate model function
                  preds <- do.call(this_model_fun,
                          args = list(params = this_coef,
                                   dose = newdata$Dose,
                                   time = newdata$Time_trans,
                                   route = newdata$Route,
                                   medium = newdata$Media))

                })
         },
         simplify = FALSE,
         USE.NAMES = TRUE)
}

#' Get residuals
#'
#' Extract residuals from a fitted `pk` object.
#'
#' Residuals are `observed - predicted`.
#'
#' @param obj A `pk` object
#' @param newdata Optional: A `data.frame` with new data for which to make
#'   predictions and compute residuals. If NULL (the default), then residuals
#'   will be computed for the data in `obj$data`. `newdata` is required to
#'   contain at least the following variables: `Time`, `Dose`, `Route`, and
#'   `Media`. `Time` will be transformed according to the transformation in
#'   `obj$scales$time` before residuals are calculated.
#' @param model Optional: Specify one or more of the fitted models for which to
#'   make predictions and calculate residuals. If NULL (the default), residuals
#'   will be returned for all of the models in `obj$stat_model`.
#' @param method Optional: Specify one or more of the [optimx::optimx()] methods
#'   for which to make predictions and calculate residuals. If NULL (the
#'   default), residuals will be returned for all of the models in
#'   `obj$optimx_settings$method`.
#' @param type Either `"conc"` (the default) or `"auc"`. `type = "conc"`
#'   predicts concentrations; `type = "auc"` predicts area under the
#'   concentration-time curve (AUC). Currently, only `type = "conc"` is
#'   implemented.
#' @return A named list of numeric matrices. There is one list element named for
#'   each model in `obj`'s [stat_model()] element, i.e. each PK model that was
#'   fitted to the data. Each list element is a matrix with the same number of
#'   rows as the data in `obj$data` (corresponding to the rows in `obj$data`),
#'   and as many columns as there were [optimx::optimx()] methods (specified in
#'   [settings_optimx()]). The column names are the method names.  Each column
#'   contains the residuals (observed - predicted) of the model fitted by the
#'   corresponding method. These residuals are concentrations in the same units
#'   as `obj$data$Conc.Units`; any concentration transformations (in
#'   `obj$scale$conc`) are *not* applied.
#' @export
#' @author Caroline Ring
residuals.pk <- function(obj,
                         newdata = NULL,
                         model = NULL,
                         method = NULL,
                         type = "conc"){
  if(!(type %in% "conc")) stop(paste("Error in residuals.pk():",
                                     "only type = 'conc' is currently implemented;",
                                     "residuals for type = 'auc' are not available yet"))
  if(is.null(model)) model <- names(obj$stat_model)
  if(is.null(method)) method <- obj$optimx_settings$method
  if(is.null(newdata)) newdata <- obj$data

  preds <- predict(obj,
                   newdata = newdata,
                   model = model,
                   method = method,
                   type = type)

  obs <- ifelse(type %in% "conc",
                obj$data$Conc,
                NA_real_)

  sapply(preds,
         function(this_pred){
            obs - this_pred
         },
         simplify = FALSE,
         USE.NAMES = TRUE)
}

#' Root mean squared error
#'
#' Extract root mean squared error of a fitted `pk` object
#'
#' @param obj A `pk` object
#' @param newdata Optional: A `data.frame` with new data for which to make
#'   predictions and compute RMSsE. If NULL (the default), then RMSEs will be
#'   computed for the data in `obj$data`. `newdata` is required to contain at
#'   least the following variables: `Time`, `Dose`, `Route`, and `Media`. `Time`
#'   will be transformed according to the transformation in `obj$scales$time`
#'   before RMSEs are calculated.
#' @param model Optional: Specify one or more of the fitted models for which to
#'   make predictions and calculate RMSEs. If NULL (the default), RMSEs will be
#'   returned for all of the models in `obj$stat_model`.
#' @param method Optional: Specify one or more of the [optimx::optimx()] methods
#'   for which to make predictions and calculate RMSEs. If NULL (the default),
#'   RMSEs will be returned for all of the models in
#'   `obj$optimx_settings$method`.
#' @param type Either `"conc"` (the default) or `"auc"`. `type = "conc"`
#'   predicts concentrations; `type = "auc"` predicts area under the
#'   concentration-time curve (AUC). Currently, only `type = "conc"` is
#'   implemented.
#' @return A named list of numeric vectors. There is one list element named for
#'   each model in `obj`'s [stat_model()] element, i.e. each PK model that was
#'   fitted to the data. Each list element is a numeric vector with as many
#'   elements as there were [optimx::optimx()] methods (specified in
#'   [settings_optimx()]). The vector names are the method names.  Each vector
#'   element contains the root mean squared error (`sqrt(mean(observed -
#'   predicted)^2)`) of the model fitted by the corresponding method, using the
#'   data in `newdata`. These RMSEs are concentrations in the same units as
#'   `obj$data$Conc.Units`; any concentration transformations (in
#'   `obj$scale$conc`) are *not* applied.
#' @export
#' @author Caroline Ring
rmse.pk <- function(obj,
                    newdata = NULL,
                    model = NULL,
                    method = NULL,
                    type = "conc"){
  if(!(type %in% "conc")) stop(paste("Error in residuals.pk():",
                                     "only type = 'conc' is currently implemented;",
                                     "residuals for type = 'auc' are not available yet"))
  if(is.null(model)) model <- names(obj$stat_model)
  if(is.null(method)) method <- obj$optimx_settings$method
  if(is.null(newdata)) newdata <- obj$data

resids <- residuals(obj,
                    newdata = newdata,
                    model = model,
                    method = method,
                    type = type)

sapply(resids,
       function(this_resid){
         apply(this_resid,
               2,
               function(x) sqrt(mean(x*x)))
       },
       simplify = FALSE,
       USE.NAMES = TRUE)
}

#' Log-likelihood
#'
#' Extract log-likelihood(s) from a fitted `pk` object
#'
#' For details on how the log-likelihood is calculated, see [log_likelihood()].
#'
#' @param obj A `pk` object
#' @param newdata Optional: A `data.frame` with new data for which to compute
#'   log-likelihood. If NULL (the default), then log-likelihoods will be
#'   computed for the data in `obj$data`. `newdata` is required to contain at
#'   least the following variables: `Time`, `Dose`, `Route`,`Media`, `Conc`,
#'   `Detect`, `Conc_SD`, `N_Subjects`. Before log-likelihood is calculated, `Time` will be
#'   transformed according to the transformation in `obj$scales$time` and `Conc`
#'   will be transformed according to the transformation in `obj$scales$conc`.
#' @param model Optional: Specify one or more of the fitted models for which to
#'   calculate log-likelihood. If NULL (the default), log-likelihoods will be
#'   returned for all of the models in `obj$stat_model`.
#' @param method Optional: Specify one or more of the [optimx::optimx()] methods
#'   for which to calculate log-likelihoods. If NULL (the default),
#'   log-likelihoods will be returned for all of the models in
#'   `obj$optimx_settings$method`.
#' @return A named list of numeric vectors. There is one list element named for
#'   each model in `obj`'s [stat_model()] element, i.e. each PK model that was
#'   fitted to the data. Each list element is a numeric vector with as many
#'   elements as there were [optimx::optimx()] methods (specified in
#'   [settings_optimx()]). The vector names are the method names.  Each vector
#'   element contains the log-likelihood of the model fitted by the
#'   corresponding method, calculated for the data in `newdata` if any has been
#'   specified, or for the data in `obj$data` if `newdata` is `NULL`.
#' @export
#' @importFrom stats logLik
#' @author Caroline Ring
logLik.pk <- function(obj,
                      newdata = NULL,
                      model = NULL,
                      method = NULL,
                      negative = FALSE,
                      force_finite = FALSE){
  if(is.null(model)) model <- names(obj$stat_model)
  if(is.null(method)) method <- obj$optimx_settings$method
  if(is.null(newdata)) newdata <- obj$data

    #transform time needed
    #first, default to identity transformation if none is specified
    if(is.null(obj$scales$time$new_units)){
      obj$scales$time$new_units <- "identity"
    }

    from_units <- unique(newdata$Time.Units)
    to_units <- ifelse(obj$scales$time$new_units %in% "identity",
                       from_units,
                       obj$scales$time$new_units)


    if(obj$scales$time$new_units %in% "auto"){
      to_units <- auto_units(y = newdata$Time,
                             from = from_units)
    }

    newdata$Time_trans <- tryCatch(convert_time(x = newdata$Time,
                                                from = from_units,
                                                to = to_units,
                                                inverse = FALSE),
                                   error = function(err){
                                     warning(paste("invivopkfit::logLik.pk():",
                                                   "Error in transforming time in `newdata` using convert_time():",
                                                   err$message))
                                     return(NA_real_)
                                   })

  #Apply concentration transformation
    newdata$Conc_trans <- rlang::eval_tidy(obj$scales$conc$expr,
                                           data = cbind(newdata,
                                                        data.frame(".conc" = newdata$Conc)))
    newdata$Conc_SD_trans <- rlang::eval_tidy(obj$scales$conc$expr,
                                           data = cbind(newdata,
                                                        data.frame(".conc" = newdata$Conc_SD)))

  sapply(model, function(this_model){
      #if there is newdata, then we have to evaluate the log-likelihood
      #get model parameters
      coefs <- coef(obj = obj,
                    model = this_model,
                    method = method)[[1]] #we have to put [[1]] because for one model,coef.pk() returns a one-element list
      #for each row of model parameters, evaluate log-likelihood
      ll <- apply(coefs,
                  1,
                  function(this_coef_row){
                    log_likelihood(par = as.list(this_coef_row),
                                   fitdata = newdata,
                                   data_sigma_group = obj$stat_error_model$data_sigma_group,
                                   modelfun = obj$stat_model[[this_model]]$conc_fun,
                                   scales_conc = obj$scales$conc,
                                   negative = negative,
                                   force_finite = force_finite)
                  })
      #set attribute "df", the number of parameters optimized for this model
      attr(ll, which = "df") <- attr(obj$stat_model[[this_model]]$fit, "npar")
      #set attributes "nobs", the number of observations in `newdata`
      attr(ll, which = "nobs") <- nrow(newdata)
      ll
    },
    simplify = FALSE,
    USE.NAMES = TRUE)

}

#' Akaike information criterion
#'
#' Get the Akaike information criterion (AIC) for a fitted `pk` object
#'
#' The AIC is calculated from the log-likelihood (LL) as follows:
#' \deqn{\textrm{AIC} = -2\textrm{LL} + k n_{par}}
#'
#' where \eqn{n_{par}} is the number of parameters in the fitted model, and
#' \eqn{k = 2} for the standard AIC.
#'
#' @param obj A `pk` object
#' @param newdata Optional: A `data.frame` with new data for which to compute
#'   log-likelihood. If NULL (the default), then log-likelihoods will be
#'   computed for the data in `obj$data`. `newdata` is required to contain at
#'   least the following variables: `Time`, `Dose`, `Route`,`Media`, `Conc`,
#'   `Detect`, `N_Subjects`. Before log-likelihood is calculated, `Time` will be
#'   transformed according to the transformation in `obj$scales$time` and `Conc`
#'   will be transformed according to the transformation in `obj$scales$conc`.
#' @param model Optional: Specify one or more of the fitted models for which to
#'   calculate log-likelihood. If NULL (the default), log-likelihoods will be
#'   returned for all of the models in `obj$stat_model`.
#' @param method Optional: Specify one or more of the [optimx::optimx()] methods
#'   for which to make predictions and calculate AICs. If NULL (the default),
#'   log-likelihoods will be returned for all of the models in
#'   `obj$optimx_settings$method`.
#' @return A named list of numeric vectors. There is one list element
#'   named for each model in `obj`'s [stat_model()] element, i.e. each PK model
#'   that was fitted to the data. Each list element is a numeric vector with as
#'   many elements as there were [optimx::optimx()] methods (specified in
#'   [settings_optimx()]). The vector names are the method names.  Each vector
#'   element contains the AIC of the model fitted by the corresponding method,
#'   using the data in `newdata`.
#' @seealso [BIC.pk(0)], [logLik.pk()]
#' @importFrom stats AIC
#' @export
#' @author Caroline Ring
AIC.pk <- function(obj,
                   newdata = NULL,
                   model = NULL,
                   method = NULL,
                   k = 2){
  if(is.null(model)) model <- names(obj$stat_model)
  if(is.null(method)) method <- obj$optimx_settings$method
  if(is.null(newdata)) newdata <- obj$data
 #get log-likliehoods
  ll <- logLik(obj = obj,
               newdata = newdata,
               model = model,
               method = method,
               negative = FALSE,
               force_finite = FALSE)
  #get number of parameters (excluding any constant, non-optimized parameters)
 AIC <- sapply(model,
                function(this_model){
                  npar <- attr(ll[[this_model]], "df")
                  -2*ll[[this_model]] + k* npar
                },
                simplify = FALSE,
                USE.NAMES = TRUE)
return(AIC)
}

#' Bayesian information criterion
#'
#' Get the Bayesian information criterion (AIC) for a fitted `pk` object
#'
#' The BIC is calculated from the log-likelihood (LL) as follows:
#' \deqn{\textrm{BIC} = -2\textrm{LL} + \log(n_{obs}) n_{par}}
#'
#' where \eqn{n_{par}} is the number of parameters in the fitted model.
#'
#' Note that the BIC is just the AIC with \eqn{k = \(n_{obs})}.
#'
#' @param obj A `pk` object
#' @param newdata Optional: A `data.frame` with new data for which to compute
#'   log-likelihood. If NULL (the default), then BICs will be
#'   computed for the data in `obj$data`. `newdata` is required to contain at
#'   least the following variables: `Time`, `Dose`, `Route`,`Media`, `Conc`,
#'   `Detect`, `N_Subjects`. Before log-likelihood is calculated, `Time` will be
#'   transformed according to the transformation in `obj$scales$time` and `Conc`
#'   will be transformed according to the transformation in `obj$scales$conc`.
#' @param model Optional: Specify one or more of the fitted models for which to
#'   calculate BIC. If NULL (the default), log-likelihoods will be
#'   returned for all of the models in `obj$stat_model`.
#' @param method Optional: Specify one or more of the [optimx::optimx()] methods
#'   for which to calculate BICs. If NULL (the default),
#'   log-likelihoods will be returned for all of the methods in
#'   `obj$optimx_settings$method`.
#' @return A named list of numeric vectors. There is one list element
#'   named for each model in `obj`'s [stat_model()] element, i.e. each PK model
#'   that was fitted to the data. Each list element is a numeric vector with as
#'   many elements as there were [optimx::optimx()] methods (specified in
#'   [settings_optimx()]). The vector names are the method names.  Each vector
#'   element contains the BIC of the model fitted by the corresponding method,
#'   using the data in `newdata`.
#' @seealso [AIC.pk(0)], [logLik.pk()]
#' @importFrom stats BIC
#' @export
#' @author Caroline Ring
BIC.pk <- function(obj,
                   newdata = NULL,
                   model = NULL,
                   method = NULL){
  if(is.null(model)) model <- names(obj$stat_model)
  if(is.null(method)) method <- obj$optimx_settings$method
  if(is.null(newdata)) newdata <- obj$data

 BIC <- AIC(obj = obj,
            newdata = newdata,
            model = model,
            method = method,
            k = log(nrow(newdata)))
  return(BIC)
}

#' Print summary of a `pk` object
#'
#' This summary includes summary information about the data; about any data
#' transformations applied; about the models being fitted; about the error model
#' being applied; and any fitting results, if the `pk` object has been fitted.
#' It also includes TK quantities calculated from the fitted model parameters,
#' e.g. halflife; clearance; tmax; Cmax; AUC; Css.
#'
#' @param obj A [pk] object.
#' @return A `data.frame` consisting of a summary table of fitting options and results.
#' @export
#' @author Caroline Ring
#'
summary.pk <- function(obj){
  suppress.messages <- obj$data_settings$suppress.messages
  #get model coefficients
  coefs <- coef(obj)
  #get coefficient SDs
  coef_sds <- coef_sd(obj)
  #transpose
  coefs <- lapply(coefs, t)
  coef_sds <- lapply(coef_sds, t)



  #For each model:
  outDF_list <- sapply(names(obj$stat_model),
         function(this_model){
           this_coef <- coefs[[this_model]]
           this_sd <- coef_sds[[this_model]]
    #loop over optimx methods (rownames of this_fit)
   outDF_model_list <- sapply(colnames(this_coef),
          function(this_method){
            #grab par_DF (with bounds & starting values)
            this_outDF <- obj$stat_model[[this_model]]$par_DF
            this_outDF$model <- this_model
            #for this method
            this_outDF$method <- this_method
            #pull fitted values

            this_outDF$fitted_value <- sapply(this_outDF$param_name,
                                              function(x) {
                                                ifelse(x %in% rownames(this_coef),
                                                       this_coef[x, this_method],
                                                       NA_real_)
                                              },
                                              simplify = TRUE,
                                              USE.NAMES = TRUE)

            this_outDF$fitted_sd <- sapply(this_outDF$param_name,
                                           function(x) {
                                             ifelse(x %in% rownames(this_sd),
                                                    this_sd[x, this_method],
                                                    NA_real_)
                                           },
                                           simplify = TRUE,
                                           USE.NAMES = TRUE)

            #get RMSE
            this_outDF$rmse <- rmse.pk(obj = obj,
                                    newdata = NULL,
                                    model = this_model,
                                    method = this_method,
                                    type = "conc")

            return(this_outDF)
    },
    simplify = FALSE,
    USE.NAMES = TRUE)

   outDF_model <- do.call(rbind, outDF_model_list)
   return(outDF_model)
  },
  simplify = FALSE,
  USE.NAMES = TRUE)

  outDF <- do.call(rbind, outDF_list)

  #add the error grouping
  outDF$error_group <- paste0("~",
                              paste(
    sapply(obj$stat_error_model$error_group,
                              rlang::as_label),
    collapse = "+")
  )

  #add the data transformations

  #add the rmse


return(outDF)
}


#' Check status of a `pk` object
#'
#' Check status of a `pk` object
#'
#' `pk` objects have integer statuses reflecting what stage of the analysis
#' process they are at.
#'
#' 1. Object has been initialized
#' 2. Data pre-processing complete
#' 3. Model pre-fitting complete 4
#' . Model fitting complete
#'
#' If a `pk` object of status 2 or greater has its instructions modified with
#' `+`, then its status will be reset to 1, indicating that any analysis results
#' contained in the object are now outdated and all steps of the analysis need
#' to be re-run.
#'
#' This function allows the user to check the status of a `pk` object.
#'
#' A message will be printed listing the analysis steps that have been completed
#' for this `pk` object, and the integer status will be returned.
#'
#' @param obj A `pk` object
#' @return The status of the `pk` object as an integer.
get_status.pk <- function(obj){
  steps <- c("1. Object has been initialized",
             "2. Data pre-processing complete",
             "3. Model pre-fitting complete",
             "4. Model fitting complete")
  message(paste(steps[seq(1, obj$status)],
                collapse = "/n"))
  return(obj$status)
}

#'Non-compartmental analysis
#'
#'Do non-compartmental analysis of a `pk` object.
#'
#'This function calls [PK::nca()] to calculate the following quantities:
#'
#'For intravenous (IV) bolus administration:
#'
#' * `nca.iv.AUC_tlast`: AUC (area under concentration-time curve) evaluated at the last reported time point
#' * `nca.iv.AUC_inf`: AUC evaluated at time t = \eqn{\infty}
#' * `nca.iv.AUMC_inf`: AUMC (area under the first moment curve, i.e. under the AUC vs. time curve) evaluated at time \eqn{t = \infty}
#' * `nca.iv.MRT`: MRT (mean residence time)
#' * `nca.iv.halflife`: Half-life
#' * `nca.iv.Clearance`: Clearance
#' * Vss: apparent volume of distribution at steady-state.
#'
#'For oral bolus administration:
#'
#' * `nca.oral.AUC_tlast`: AUC (area under concentration-time curve) evaluated at the last reported time point
#' * `nca.oral.AUC_inf`: AUC evaluated at time \eqn{t = \infty}
#' * `nca.oral.AUMC_inf`: AUMC (area under the first moment curve, i.e. under the AUC vs. time curve) evaluated at time \eqn{t = \infty}
#' * `nca.oral.MTT`: MTT (mean transit time), the sum of MRT and mean absorption time (MAT)
#'
#'If `obj$data` only contains data for one of those routes, the NCA quantities
#'for the "missing" route will be returned, but filled with `NA_real_.`
#'
#'Additionally, for oral bolus administration, the following quantities are
#'estimated using [get_peak()]:
#'
#' * `nca.oral.tmax`: the time at which peak concentration occurs
#' * `nca.oral.Cmax`: the peak concentration
#'
#'If `dose_norm == TRUE`, then all data are dose-normalized before computing NCA
#'quantities. This means that all dose-dependent NCA quantities reflect an
#'estimate for a 1 mg/kg dose.
#'
#' If `dose_norm == FALSE`, then all quantities are calculated separately for each dose.
#'
#' @param obj A `pk` object containing concentration-time data in the element `obj$data`.
#' @param dose_norm Whether to dose-normalize data before performing NCA. Default TRUE.
#' @return A `data.frame` with variables as listed in Details.
nca.pk <- function(obj,
                   newdata = NULL,
                   dose_norm = TRUE){

  if(is.null(newdata)){
    newdata <- obj$data
  }

  nca_names <- c("AUC_tlast",
                 "AUC_inf",
                 "AUMC_inf",
                 "MRT",
                 "halflife",
                 "Clearance",
                 "Vss")

  if("iv" %in% newdata$Route){
    iv_data <- subset(newdata,
                      Route %in% "iv")
    if(dose_norm %in% TRUE){
      nca_iv <- do_nca(obs = iv_data,
                       dose_norm = dose_norm)
    }else{
      #do NCA separately for each dose
      dose_list <- split(iv_data,
                         iv_data$Dose)
      nca_iv <- do.call(rbind,
                        sapply(dose_list,
                               function(this_data){
                                 this_nca <- do_nca(obs = this_data,
                                                    dose_norm = FALSE)
                                 this_nca <- cbind(Dose = this_data$Dose)
                                 this_nca
                               })
      )
    }
    names(nca_iv) <- paste0("nca.iv.",
                            names(nca_iv))
  }else{ #if no IV data, fill with NAs
    nca_iv <- data.frame(as.list(rep(NA_real_,
                                     length(nca_names))))
    names(nca_iv) <- paste0("nca.iv.",
                            nca_names)
  }

  if("po" %in% newdata$Route){
    oral_data <- subset(newdata,
                        Route %in% "po")
    if(dose_norm %in% TRUE){
      nca_oral <- do_nca(obs = oral_data,
                          dose_norm = dose_norm)
      nca_oral <- cbind("nca.oral.Dose" = 1,
                        nca_oral)
      max_df <- as.data.frame(
        get_peak(x = oral_data$Time,
                 y = oral_data$Conc/oral_data$Dose)
      )
      names(max_df) <- c("tmax",
                         "Cmax")
      max_df <- cbind("Dose" = 1,
                      max_df)
    }else{
      #do NCA separately for each dose
      dose_list <- split(oral_data,
                         oral_data$Dose)
      nca_oral <- do.call(rbind,
                          sapply(dose_list,
                                 function(this_data){
                                   this_nca <- do_nca(obs = this_data,
                                                      dose_norm = FALSE)
                                   this_nca <- cbind(Dose = this_data$Dose)
                                   this_nca
                                 })
      )
      max_df <- do.call(rbind,
                        sapply(dose_list,
                               function(this_data){
                                 max_list <- as.data.frame(
                                   get_peak(x = this_data$Time,
                                            y = this_data$Conc/this_data$Dose)
                                 )
                                 names(max_list) <- c("tmax",
                                                      "Cmax")
                                 max_list <- cbind("Dose" = this_data$Dose,
                                                   max_list)
                                 max_list
                               }
                        )
      )
    }
    names(nca_oral) <- paste0("nca.oral.",
                              names(nca_oral))

    nca_oral <- merge(nca_oral,
                      max_list,
                      by = "Dose")
  }else{
    nca_oral <- data.frame(as.list(rep(NA_real_,
                                       length(nca_names) + 2)))
    names(nca_oral) <- paste0("nca.oral.",
                              c("Dose",
                                nca_names,
                                c("tmax",
                                  "Cmax")
                              )
    )
  }

  names(nca_oral)[match(c("nca.oral.Clearance",
                          "nca.oral.MRT"))] <- c("nca.oral.Clearance_Fgutabs",
                                                 "nca.oral.MTT")
  #remove oral halflife and Vss estimates because they are not valid
  nca_oral[c("nca.oral.halflife",
             "nca.oral.Vss")] <- NULL

  nca_out <- cbind(nca_iv,
                   nca_oral)
  return(nca_out)
}

#' Plot concentration vs. time data.
#'
#' @param obj A [pk()] object with concentration-time data in element `obj$data`.
#' @param newdata Optional: A new set of concentration vs. time data to plot,
#'   different from `obj$data`. Default `NULL` to plot the data in `obj$data`.
#' @return A [ggplot2::ggplot()] plot object.
#' @export
#' @author Caroline Ring
plot_data.pk <- function(obj,
                         newdata = NULL,
                         plot_dose_norm = TRUE,
                         plot_log10_conc = FALSE){
  if(is.null(newdata)){
    newdata <- obj$data
  }

  #set the Detect column to a factor with levels "Detect" and "Non-Detect"
  if(!("Detect" %in% names(newdata))){
    newdata$Detect <- TRUE
  }

  newdata$Detect <- factor(newdata$Detect,
                           levels = c(TRUE, FALSE),
                           labels = c("Detect", "Non-Detect"))

  if(!("Conc_SD" %in% names(newdata))){
    newdata$Conc_SD <- 0
  }

  #create generic named columsn for plotting, containing either dose-normalized or
  #non-dose-normalized concentrations, depending on `plot_dose_norm`
  if(plot_dose_norm %in% TRUE){
    newdata$conc_plot <- newdata$Conc_Dose
    newdata$conc_sd_plot <- newdata$Conc_SD_Dose
  }else{
    newdata$conc_plot <- newdata$Conc
    newdata$conc_sd_plot <- newdata$Conc_SD
  }


  #generate plot title
  plot_title <- paste0(obj$data_info$Chemical,
                       " (", unique(obj$data$Compound), ")\n",
                       "Species = ", obj$data_info$Species, ", ",
                       "Doses = ", paste(signif(
                         sort(unique(obj$data$Dose)),
                         3),
                         collapse = ", "), " mg/kg\n")
  #now plot
  p <- ggplot(newdata,
              aes(x = Time,
                  y = conc_plot)) +
    geom_blank()

  if(plot_log10_conc %in% FALSE){ #plot error bars
    if(plot_dose_norm %in% TRUE){  #map colors to dose
      p <- p +
        geom_errorbar(aes(ymin = conc_plot - conc_sd_plot,
                          ymax = conc_plot + conc_sd_plot,
                          color = Dose))
    }else{   #do not map color to dose
      p <- p +
        geom_errorbar(aes(ymin = conc_plot - conc_sd_plot,
                          ymax = conc_plot + conc_sd_plot))
    }
  } #end if(plot_log10_conc %in% FALSE)

  if(plot_dose_norm %in% TRUE){ #mapping color to dose
    #concentration-dose observation points:
    #shape mapped to Reference, color to Dose, fill yes/no to Detect
    #first plot detected points with white fill
    p <- p +
      #plot detected points with white fill
      geom_point(data = newdata[Detect %in% "Detect"],
                 aes(shape = Reference,
                     color = Dose),
                 fill = "white",
                 size = 4,
                 stroke = 1.5) +
      #plot detected points with fill mapped to dose, but alpha mapped to Detect
      geom_point(data = newdata[Detect %in% "Detect"],
                 aes(shape = Reference,
                     color = Dose,
                     fill = Dose,
                     alpha = Detect),
                 size = 4,
                 stroke = 1.5) +
      #then plot non-detect points with white fill,
      #and position jittered
      geom_jitter(data = newdata[Detect %in% "Non-Detect"],
                  aes(shape = Reference,
                      color = Dose),
                  fill = "white",
                  size = 4,
                  stroke = 1.5,
                  width = jitter_nondetect,
                  height = 0) +
      #then plot non-detect points with fill, but alpha mapped to Detect
      #and position jittered
      geom_jitter(data = newdata[Detect %in% "Non-Detect"],
                  aes(shape = Reference,
                      color = Dose,
                      fill = Dose,
                      alpha = Detect),
                  size = 4,
                  stroke = 1.5,
                  width = jitter_nondetect, height = 0)

    p <- p +
      facet_grid(rows = vars(Route),
                 cols = vars(Media),
                 scales = "free_y")

    p <- p +
      scale_color_viridis_c(name = "Dose, mg/kg") +
      scale_fill_viridis_c(na.value = NA, name = "Dose, mg/kg")

  }else{ #if plot_dose_norm %in% FALSE, don't map color to dose
    #concentration-dose observation points:
    #shape mapped to Reference, fill yes/no to Detect
    #first plot detected points with white fill
    p <- p +
      #plot detected points with white fill
      geom_point(data = newdata[Detect %in% "Detect"],
                 aes(shape = Reference),
                 color = "gray50",
                 fill = "white",
                 size = 4,
                 stroke = 1.5) +
      #plot detected points with alpha mapped to Detect
      geom_point(data = newdata[Detect %in% "Detect"],
                 aes(shape = Reference,
                     alpha = Detect),
                 color = "gray50",
                 fill = "gray50",
                 size = 4,
                 stroke = 1.5) +
      #then plot non-detect points with white fill,
      #and position jittered
      geom_jitter(data = newdata[Detect %in% "Non-Detect"],
                  aes(shape = Reference),
                  color = "gray50",
                  fill = "white",
                  size = 4,
                  stroke = 1.5,
                  width = jitter_nondetect,
                  height = 0) +
      #then plot non-detect points with fill, but alpha mapped to Detect
      #and position jittered
      geom_jitter(data = newdata[Detect %in% "Non-Detect"],
                  aes(shape = Reference,
                      alpha = Detect),
                  color = "gray50",
                  fill = "gray50",
                  size = 4,
                  stroke = 1.5,
                  width = jitter_nondetect, height = 0)

    #do facet_wrap by combinations of route, media, dose
    p <- p + facet_wrap(vars(Route, Media, Dose),
                        labeller = "label_both",
                        scales = "free")
  } #end check if plot_dose_norm %in% TRUE/FALSE

  p <- p +
    scale_shape_manual(values = 21:25) + #use only the 5 shapes where a fill can be added
    #this limits us to visualizing only 5 References
    #but admittedly it's hard to distinguish more than 5 shapes anyway
    #if detect =FALSE, fully transparent; if detect = TRUE, fully solid
    scale_alpha_manual(values = c("Detect" = 1,
                                  "Non-Detect" = 0),
                       drop = FALSE,
                       name = NULL) +
    guides(alpha = guide_legend(override.aes = list(shape = 21,
                                                    color = "black",
                                                    stroke = 1,
                                                    fill = c("black", NA),
                                                    alpha = 1)
    )
    ) +
    labs(title = plot_title) +
    xlab("Time, hr")

  if(plot_dose_norm %in% TRUE){
    p <- p + ylab("Concentration/Dose, (mg/kg)/(mg/L)")
  }else{
    p <- p + ylab("Concentration, mg/L")
  }

  p <- p  +
    theme_bw() +
    theme(plot.title = element_text(size = 12),
          strip.background = element_blank())
  #log-scale y axis if so specified
  if(plot_log10_conc %in% TRUE){
    p <- p + scale_y_log10() + annotation_logticks(sides = "l")
  }

  return(p)

}

#' Rescale time
#'
#' If requested, rescale time from hours to some more convenient unit.
#'
#' Time is rescaled based upon the time of the last detected observation in units of hours.
#'
#'  | Range of Last Detect Time in Hours | New Time Units | Range of New Last Detect Time |
#'  | ----------- | -------------- |
#'  | \\[0, 0.5) | minutes | \\[0, 30) |
#'  | \\[0.5, 24) | hours | \\[0.5, 24) |
#'  | \\[24, 720) | days | \\[1, 30) |
#'  | \\[720, 8766) | months | \\[0.99, 12) |
#'  | \\[8766, 87660) | years | \\[1, 10) |
#'  | \\[87660, 876600) | decades | \\[1, 10) |
#'  | \\[876600, Inf) | centuries | \\[1, Inf) |
#'
#'  (At present, no data in CvT requires decades or centuries)
#
#' @param obj A [pk()] object.
#' @return The `data` element of the object, but modified: the original
#'   time (in hours) is now in variable `Time.Hours`; rescaled time is in
#'   `Time`; units of rescaled time are in `Time.Units`.
#' @author Caroline Ring
rescale_time.pk <- function(obj){
  #Rescale time if so requested
  #Save the original time in units of hours
  obj$data$Time.Hours <- obj$data$Time
  #Now do the rescale
  if(obj$data_trans$rescale_time %in% TRUE){
    last_detect_time <- obj$data_info$last_detect_time

    new_time_units <- cut(last_detect_time,
                          breaks = c(0,
                                     0.5,
                                     24,
                                     720,
                                     8766,
                                     87660,
                                     876600,
                                     Inf),
                          labels = c("minutes",
                                     "hours",
                                     "days",
                                     "months",
                                     "years",
                                     "decades",
                                     "centuries"),
                          right = FALSE)

    obj$data$Time <- convert_time(x = obj$data$Time.Hours,
                                  from = "hours",
                                  to = new_time_units,
                                  inverse = FALSE)

  }else{
    new_time_units <- "hours"
  }

  obj$data$Time.Units <- new_time_units

  return(obj$data)
}



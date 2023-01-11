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
#' | sigma          | 1e-4         | Default         |
#'
#' If `Vdist_from_species == TRUE`, then this function will set a lower bound
#' for `Vdist` or `V1` equal to 1/10 of the physiological plasma volume for the
#' species represented in `fitdata`. The physiological plasma volume is taken
#' from [httk::physiology.data], and represents the average plasma volume (mL/kg
#' body weight) multiplied by the average body weight (kg), converted from mL to
#' L by dividing by 1000.
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
#'@param  Vdist_from_species Logical: TRUE to estimate the lower bound
#'   of `Vdist` or `V1` as 1/10 of species-specific physiological plasma
#'   volume. FALSE to use the lower bound for `Vdist` or `V1`
#'   specified in `lower_default`.
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
                                               1e-4 #sigma
                                               ),
                               lower_bound_msg = "Default"
                             ),
                             Vdist_from_species = FALSE,
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

  if(Vdist_from_species %in% TRUE &
     model %in% c("1compartment", "2compartment")){
    #For Vdist or V1: Set the theoretical lower bound to something on the order of
    #the species-specific total plasma volume, pulled from httk physiology.data
    phys <- httk::physiology.data
    species <- unique(tolower(fitdata$Species)) #there should be only one species
    #if httk::physiology.data has data on this species:
    if(any(tolower(names(phys)) %in%  species)){
      #get the column number corresponding to this species
      species_col <- which(tolower(names(phys)) %in%
                             species)
      #pull the average plasma volume in mL/kg
      plasma_vol_mL_kg <- phys[phys$Parameter %in% "Plasma Volume",
                               species_col]
      #pull the average body weight in kg
      body_wt_kg <- phys[phys$Parameter %in% "Average BW",
                         species_col]
      #convert mL/kg to L plasma
      plasma_vol_L <- plasma_vol_mL_kg * body_wt_kg * 1e-3
      #assuming this gives us a valid answer:
      if(is.finite(plasma_vol_L)){
        #assign 1/10 of this value as theoretical lower bound on Vdist
        par_DF[par_DF$param_name %in%
                 c("Vdist", "V1"),
               "lower_bound"] <- plasma_vol_L/10
        #add message explaining lower bound & source
        par_DF[par_DF$param_name %in%
                 c("Vdist", "V1"),
               "lower_bound_msg"] <- paste0("0.1 * ",
                                            "average total plasma volume for species (",
                                            signif(plasma_vol_L, digits = 3),
                                            " L), from httk::physiology.data")
      } #end if(is.finite(plasma_vol_L))
    } #end if(any(tolower(names(phys)) %in%  species))
  } #end if if(Vdist_from_species %in% TRUE)

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

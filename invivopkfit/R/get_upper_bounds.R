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
#' | Fgutabs_Vdist  | 1e4    | Default         |
#' | Fgutabs_V1     | 1e4    | Default         |
#' | sigma          | Inf         | Default         |
#'
#' #' Any parameters which will not be estimated from data (based on either the
#' variable `optimize_param` in `par_DF` if `par_DF` is provided, or the output
#' of [get_opt_params()] if `par_DF` is not provided) are assigned an upper bound
#' of `NA_real_`.
#'
#' @param fitdata A `data.frame`: the concentration-time-dose data to be used for
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
#' @param Fgutabs_Vdist_from_species Logical: TRUE to estimate the upper bound
#'   of `Fgutabs_Vdist` or `Fgutabs_V1` using the upper bound of `Fgutabs`
#'   specified in `upper_default`, divided by a species-specific lower bound on
#'   the volume of distribution, based on species-specific physiological plasma
#'   volumes. FALSE to use the upper bound for `Fgutabs_Vdist` or `Fgutabs_V1`
#'   specified in `upper_default`.
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
                                               1e4, #Fgutabs_Vdist
                                               1e4, #Fgutabs_V1
                                               Inf), #sigma
                               upper_bound_msg = "Default"
                             ),
                             Fgutabs_Vdist_from_species = TRUE,
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


    #assign sigmas
    par_DF[grepl(x = par_DF$param_name,
                 pattern = "sigma"),
           "upper_bound"] <- upper_default[upper_default$param_name %in% "sigma",
                                           "upper_bound"]
    par_DF[grepl(x = par_DF$param_name,
                 pattern = "sigma"),
           "upper_bound_msg"] <- upper_default[upper_default$param_name %in% "sigma",
                                               "upper_bound_msg"]

    if(Fgutabs_Vdist_from_species %in% TRUE &
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
        Vdist_lower <- plasma_vol_L/10
        Fgutabs_upper <- par_DF[par_DF$param_name %in% "Fgutabs",
                                "upper_bound"]
        #assuming this gives us a valid answer:
        if(is.finite(Vdist_lower)){
          #assign 1/10 of this value as theoretical lower bound on Vdist
          par_DF[par_DF$param_name %in%
                   c("Fgutabs_Vdist", "Fgutabs_V1"),
                 "upper_bound"] <- Fgutabs_upper/Vdist_lower
          #add message explaining lower bound & source
          par_DF[par_DF$param_name %in%
                   c("Fgutabs_Vdist", "Fgutabs_V1"),
                 "upper_bound_msg"] <- paste0("[Upper bound of Fgutabs,",
                                              signif(Fgutabs_upper, digits = 3),
                 "]/ [lower bound of Vdist, 0.1 * ",
                                              "average total plasma volume for species (",
                                              signif(plasma_vol_L, digits = 3),
                                              " L), from httk::physiology.data]")
        } #end if(is.finite(plasma_vol_L))
      } #end if(any(tolower(names(phys)) %in%  species))
    } #end if if(Fgutabs_Vdist_from_species %in% TRUE)

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

#' Get TK stats
#'
#' Extract derived TK statistics from a fitted [pk()] model object.
#'
#' After fitting model parameters (e.g. elimination rate, volume of
#' distribution, absorption rate, bioavailability), it can be useful to derive
#' summary toxicokinetic statistics such as total clearance rate, half-life,
#' peak concentration, AUC_inf (the area under the concentration-time curve when
#' time goes to infinity), etc.
#'
#' Many of these TK statistics depend not only on chemical and species, but also
#' on route, media (tissue), and dose. Therefore, TK stats need to be computed
#' for a specific set of Chemical, Species, Route, Media, and Dose.
#'
#' TK statistics for a defined [pk_model()] object are computed using the
#' function named in the model's `tkstats_fun`. For the built-in models, the
#' `tkstats_fun` functions are the following. See the documentation for the
#' individual functions for details on what TK stats are calculated for each
#' model, and how they are calculated.

#' -`model_1comp`: [tkstats_1comp()]
#' -`model_2comp`: [tkstats_2comp()]
#' -`model_flat`: [tkstats_flat()]
#'
#' @param obj A [pk()] model object. Must be fitted, or the function will exit
#'   with an error.
#' @param newdata Optional: A `data.frame` containing new data for which to
#'   compute the TK stats. Must contain at least variables `Chemical`,
#'   `Species`, `Route`, `Media`, `Dose`, `Dose.Unit`, `Conc.Unit`, either
#'   `Time_trans.Unit` or `Time.Unit`, and any other variables named in
#'   `tk_grouping`. Default `NULL`, to use the data in `obj$data`.
#' @param tk_group A list of variables provided using a `dplyr::vars()` call.
#'   The data (either `newdata` or `obj$data`) will be grouped according to the
#'   unique combinations of these variables. For each unique combination of
#'   these variables in the data, a set of TK statistics will be computed. The
#'   default is `obj$settings_data_info$nca_group`, to derive TK statistics for
#'   the same groups of data as non-compartmental analysis statistics. With the
#'   default, you can directly compare e.g. a model-predicted AUC_inf to the
#'   corresponding NCA-estimated AUC_inf. However, you may specify a different
#'   data grouping if you wish. Each group should have a unique combination of
#'   `Chemical`, `Species`, `Route`, `Media`, and `Dose`.
#' @param model Character: One or more of the models fitted. Default `NULL` to
#'   return TK stats for all models.
#' @param method Character: One or more of the [optimx::optimx()] methods used.
#'   Default `NULL` to return TK stats for all methods.
#' @return A `list` of `data.frame` objects, one  named for each model in
#'   `model`. Each `data.frame` will have the variables in the `data.frame`
#'   returned by the `tkstats_fun` for its corresponding model. (For the
#'   built-in models `model_flat`, `model_1comp`, and `model_2comp`, these
#'   variables are `param_name` and `param_value`.) Additionally, there will be
#'   a variable `method` denoting the [optimx::optimx()] method used to optimize
#'   the set of model parameters used to derive each set of TK statistics.
#' @export
#' @author Caroline Ring
get_tkstats.pk <- function(obj,
                           newdata = NULL,
                           tk_group = obj$settings_data_info$nca_group,
                           model = NULL,
                           method = NULL,
                           vol_unit = "L"){

  #ensure that the model has been fitted
  check <- check_required_status(obj = obj,
                                 required_status = status_fit)
  if(!(check %in% TRUE)){
    stop(attr(check, "msg"))
  }

  if(is.null(model)) model <- names(obj$stat_model)
  if(is.null(method)) method <- obj$settings_optimx$method
  if(is.null(newdata)) newdata <- obj$data

  #check that all methods are valid
  if(!(all(method %in% obj$settings_optimx$method))){
    stop(paste("All values in `method` must be found in `obj$settings_optimx$method.",
               paste0("`method` = ", paste(method, sep = ", ")),
               paste0("`obj$settings_optimx$method` = ", paste(obj$settings_optimx$method)),
               sep = "\n"))

  }

  #check that tk_group is valid: it must produce groups with a unique
  #combination of Chemical, Species, Route, Media, and Dose

  newdata_grouped <- do.call(dplyr::group_by,
                             c(list(newdata),
                               tk_group)) %>%
    dplyr::distinct(Chemical, #for each tk_group, take distinct rows by these variables
                    Species,
                    Route,
                    Media,
                    Dose) %>%
    dplyr::count() #how many distinct rows per group?

  #if more than one distinct row per group, stop
  if(any(newdata_grouped$n > 1)){
   stop("tk_group does not produce groups with unique combinations of Chemical, Species, Route, Media, and Dose.")
  }

  all_coefs <- coef(obj,
                     model = model,
                     method = method)

tkstats_all <- sapply(model,
       function(this_model){
         #get the model's TKstats function
         this_tkstats_fun <- obj$stat_model[[this_model]]$tkstats_fun
         #get any additional arguments to the model's TKstats function
         this_tkstats_args <- obj$stat_model[[this_model]]$tkstats_fun_args
         #Get the matrix of coefficients for this model -- one row named for each method
         this_coef <- all_coefs[[this_model]]
         #Derive dose units
         dose_unit <- unique(newdata$Dose.Units)
         #Derive conc unit
           conc_unit <- unique(newdata$Conc.Units)
         #Derive time unit
         if("Time_trans.Units" %in% names(newdata)){
           time_unit <- unique(newdata$Time_trans.Units)
         }else{
           time_unit <- unique(newdata$Time.Units)
         }
         #loop over methods: get tkstats for the set of coefficients for each method
         #the result will be a named list of data.frames, one for each method
         tkstats_list <- sapply(method,
                function(this_method){
                  #pull the set of coefficients for this method
                  #convert from one row of a data.frame to a named vector
                  coef_row <- unlist(this_coef[this_method, ])
                  #get the unique combinations of refrence, route, media, dose from obj$data_info$nca
                  tkstats_this_method <- do.call(dplyr::group_by, #use do.call() because the grouping is a *list* of quosures
                          c(list(newdata),
                            tk_group
                          )
                  ) %>%
                    dplyr::summarise( #for each group: get the data.frame of tkstats
                      do.call(this_tkstats_fun,
                              args = c(
                                list(pars = coef_row, #a named numeric vector
                                     route = unique(Route),
                                     medium = unique(Media),
                                     dose = unique(Dose),
                                     time_unit = time_unit,
                                     conc_unit = conc_unit,
                                     vol_unit = vol_unit
                                ),
                                this_tkstats_args
                              )
                      )#end do.call(this_tkstats_fun)
                      ) %>%
                        as.data.frame() #convert from tibble back to data.frame

                  #add method as a variable
                  tkstats_this_method$method <- this_method
                  #return the data.frame of tkstats for this model & this method
                  return(tkstats_this_method)
                },
                simplify = FALSE, #return
                USE.NAMES = TRUE #name the list elements after the items in "method"
                )
        #now, rowbind the tkstats data.frames for each method for this model
         #the result will be one data.frame
         tkstats_this_model <- do.call(rbind,
                 tkstats_list)
         rownames(tkstats_this_model) <- NULL
         #return the tkstats data.frame for this model
         return(tkstats_this_model)
       }, #end function(this_model)
       simplify = FALSE, #return a list
       USE.NAMES = TRUE #name the list elements after the items in "model"
) #end sapply over models
#the result will be a named list of data.frames, one for each model
#return it

return(tkstats_all)
}

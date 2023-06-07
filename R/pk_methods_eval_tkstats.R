eval_tkstats <- function(obj,
                         newdata = NULL,
                         model = NULL,
                         method = NULL,
                         tk_group = NULL,
                         exclude = TRUE,
                         dose_norm = TRUE,
                         which_params = NULL,
                         ...){

  #ensure that the model has been fitted
  check <- check_required_status(obj = obj,
                                 required_status = status_fit)
  if(!(check %in% TRUE)){
    stop(attr(check, "msg"))
  }

  if(is.null(model)) model <- names(obj$stat_model)
  if(is.null(method)) method <- obj$settings_optimx$method
  if(is.null(newdata)) newdata <- obj$data
  if(is.null(tk_group)) tk_group <- obj$settings_data_info$nca_group

  method_ok <- check_method(obj = obj, method = method)
  model_ok <- check_model(obj = obj, model = model)


  grp_vars <- sapply(tk_group,
                     rlang::as_label)

  newdata_ok <- check_newdata(newdata = newdata,
                              olddata = obj$data,
                              req_vars = union(
                                c("Chemical",
                                  "Species",
                                  "Time",
                                  "Time.Units",
                                  "Dose",
                                  "Conc",
                                  "Dose.Units",
                                  "Conc.Units",
                                  "Route",
                                  "Media"),
                                grp_vars),
                              exclude = exclude)

  #if exclude = TRUE, remove excluded observations
  if(exclude %in% TRUE){
    newdata <- subset(newdata, exclude %in% FALSE)
  }

#calc NCA for newdata
   nca_df <- get_nca(obj = obj,
                        newdata = newdata,
                        nca_group = tk_group,
                        exclude = exclude,
                        dose_norm = dose_norm)

   nca_df <- nca_df %>%
      tidyr::pivot_longer(cols = !(tidyselect::all_of(grp_vars)),
                                    names_to = "param_name",
                                    values_to = "param_value")

    #get tkstats
    if(dose_norm %in% TRUE){
    newdata$Conc <- newdata$Conc/newdata$Dose
    newdata$Dose <- newdata$Dose/newdata$Dose
    }

    tkstats_list <- get_tkstats(obj = obj,
                                newdata = newdata,
                                model = model,
                                method = method,
                                tk_group = tk_group,
                                exclude = exclude)

    #merge
    tk_eval <- sapply(tkstats_list,
                      function(this_tkstats){

                        tmp <- dplyr::inner_join(this_tkstats,
                                nca_df,
                                 by = c(grp_vars,
                                        "param_name"),
                                suffix = c(".fitted", ".nca")) %>%
                          dplyr::select(tidyselect::all_of(grp_vars),
                                        method,
                                        param_name,
                                        param_units,
                                        param_value.nca,
                                        param_value.fitted)

                        if(!is.null(which_params)){
                          tmp <- tmp %>%
                            dplyr::filter(param_name %in% which_params)
                        }
                        tmp
                      },
                      simplify = FALSE,
                      USE.NAMES = TRUE)


    tk_eval


}

#' Create a `pkanalysis` object from three `pkfit` objects
#'
#' Take three `pkfit` objects that differ only in which model was fit, and
#' produce a combined `pkanalysis` object that allows fit comparison across
#' models.
#'
#' @param obj_flat A `pkfit` object for the flat model
#' @param obj_1comp A `pkfit` object for the 1-compartment model
#' @param obj_2comp A `pkfit` object for the 2-compartment model
pkanalysis <- function(obj_flat,
                       obj_1comp,
                       obj_2comp){
  #Check whether all three are valid pkfit objects
  if(!(is.pkfit(obj_flat) &
    is.pkfit(obj_1comp) &
    is.pkfit(obj_2comp))){
    stop("One or more of the inputs is not a valid pkfit object.")
  }
  #Check to make sure everything except model is the same for all three objects
  expected_names <- c("data",
                      "data_info",
                      "data_trans",
                      "analysis_type",
                      "model_info",
                      "fit_params",
                      "fit_gof",
                      "fit_args",
                      "fit_control",
                      "fit_info")
  check <- sapply(setdiff(expected_names,
                 c("model_info",
                   "fit_params",
                   "fit_gof",
                   "fit_info")),
         function(this_item){
           all(all.equal(obj_flat[[this_item]],
                     obj_1comp[[this_item]]) %in% TRUE) &
             all(all.equal(obj_flat[[this_item]],
                       obj_2comp[[this_item]]) %in% TRUE)
         },
         simplify = TRUE,
         USE.NAMES = TRUE)

  if(!(all(check %in% TRUE))){
    stop(paste("The three input pkfit objects have differences in the following elements:",
         paste(names(check)[!(check %in% TRUE)],
               sep = ", ")))
  }

 #Continue with combininng
  obj <- list(
    "flat" = obj_flat,
    "1compartment" = obj_1comp,
    "2compartment" = obj_2comp
  )

  class(obj) <- c(class(obj), 'pkanalysis')

  return(obj)
}

#' Get winning model
#'
#' Get the winning model (lowest AIC)
#'
#' The "winning" model for a given analysis is the one with the lowest Akaike
#' Information Criterion (AIC).
#' @param obj
#' @return Character: The name of the winning model ("flat", "1compartment",
#'   "2compartment").
#' @export
#' @author Caroline Ring
winmodel.pkanalysis <- function(obj){
  AICs <- sapply(obj,
                 function(this_pkfit){
                   AIC.pkfit(obj= this_pkfit)
                 },
                 simplify = TRUE,
                 USE.NAMES = TRUE)
  if(!(all(is.na(AICs)))){
  winmodel <- names(AICs)[which.min(AICs)]
  }else{
    winmodel <- "No fit"
  }

  return(winmodel)
}

#' Relative likelihood of each model
#'
#' Calculate relative likelihood of the flat model vs. each model
#'
#' Relative likelihood of model 1 vs. model 2 (fitted to the same data with the same
#' transformations) is given by
#'
#' \deqn{\textrm{RL} = \exp(
#' \frac{\textrm{AIC}_2 - \textrm{AIC}_1)}{2}
#' )}
#'
#' @param obj A `pkanalysis` object
#' @return A named numeric vector: names are "flat", "1compartment", and
#'   "2compartment", each containing the relative likelihood of the flat model
#'   compared to the corresponding named model.
#' @export
#' @author Caroline Ring
rel_like.pkanalysis <- function(obj){
  AICs <- sapply(obj,
                 function(this_pkfit){
                   AIC.pkfit(obj = this_pkfit)
                 },
                 simplify = TRUE,
                 USE.NAMES = TRUE)
  rel_1comp <- exp((AICs["1compartment"] - AICs["flat"])/2)
  rel_2comp <- exp((AICs["2compartment"] - AICs["flat"])/2)
  rel_flat <- 1

  return(data.frame(model = c("flat",
                               "1compartment",
                               "2compartment"),
                     rel_like = c(rel_flat,
                                  rel_1comp,
                                  rel_2comp)))
}

plot.pkanalysis <- function(obj,
                            newdata = NULL,
                            dose_norm = NULL,
                            log_trans = NULL,
                            n_interp_time= 10){
  #get the predictions by each model
  pred_DF <- do.call(rbind,
                          sapply(obj,
                         function(this_pkfit){
                           this_pred <- get_predDF(obj = this_pkfit,
                                      newdata = newdata,
                                      dose_norm = dose_norm,
                                      n_interp_time = n_interp_time)
                           this_pred <- cbind("model" = this_pkfit$model_info$model)
                           this_pred
                         },
                         simplify = FALSE,
                         USE.NAMES = TRUE)
  )

  #plot concentration vs. time data
  p <- plot_data.pkfit(obj = obj,
                       newdata = newdata,
                       dose_norm = dose_norm,
                       log_trans = log_trans)

  #Add model fit curves
  #Map linetype to model
  p <- p +
    geom_line(data = pred_DF,
                     aes(linetype = model))

  #Update plot title/subtitle to reflect *winning* model

  #first get the name of the winning model
  winning <- winmodel.pkanalysis(obj)

  #then get the parameters of the winning model
  par <- coef.pkfit(obj[[winning]])
  #paste into a comma-separated list
  par_char <- paste(
    paste(names(par),
          signif(par, 3), #keep 3 sigfigs
          sep=" = "),
    collapse = ", ")

  plot_subtitle <- paste0("Winning Model = ",
                          winmodel)
  plot_subtitle <- paste0(plot_subtitle, "\n",
                          par_char)

  #generate plot title
  plot_title <- paste0(obj$flat$data_info$DTXSID,
                       " (", unique(obj$flat$data$Compound), ")\n",
                       "Species = ", obj$flat$data_info$Species, ", ",
                       "Doses = ", paste(signif(
                         sort(unique(obj$flat$data$Dose)),
                         3),
                         collapse = ", "), " mg/kg\n",
                       "Analysis Type = ", obj$flat$analysis_type, "\n",
                       "Fitting options: ",
                       "log-transform ",
                       obj$flat$data_trans$fit_log_conc,
                       "; dose-normalize ",
                       obj$flat$data_trans$fit_conc_dose,
                       "; rescale time ",
                       obj$flat$data_trans$rescale_time)

    #replace the existing title with the new title/subtitle
    p <- p + labs(title = plot_title,
         subtitle = plot_subtitle)
#return the ggplot2 object
    return(p)
}

#' Get summary of pkanalysis object
#'
#' @param obj
#' @return A named list
summary.pkanalysis <- function(obj,
                               dose_norm = NULL,
                               log_trans = NULL,
                               match_nondetect = TRUE){
summary_list <- sapply(obj,
       function(this_pkfit){
         summary.pkfit(this_pkfit,
                       dose_norm = dose_norm,
                       log_trans = log_trans,
                       match_nondetect = match_nondetect)
       },
       simplify = FALSE,
       USE.NAMES = TRUE
      )

#rowbind each element of the elements of summary_list
#they all have the same variables except for params, params_sd, and tk_stats
#those will have different names for each model

#for those ones, reshape longer, with a "variable" and "value" column
summary_element_names <- unique(sapply(summary_list,
                                  names))

summary_combined <- sapply(summary_element_names,
       function(this_name){
         do.call(rbind,
                 lapply(names(summary_list),
                        function(this_model){
                          this_item <- summary_list[[this_model]][[this_name]]
                          if(this_name %in% c("params",
                                              "params_sd",
                                              "tk_stats")){
                            #reshape longer
                            this_item <- data.frame("variable" = names(this_item),
                                                    "value" = unlist(this_item[1, ])
                            )
                          }

                          this_item <- cbind(
                            model = this_model,
                            this_item
                          )

                          return(this_item)
                        }
                        )
                 )
       },
       simplify = FALSE,
       USE.NAMES = TRUE
       )

#Also add the winning model and relative likelihood to the combined summary
summary_combined <- c(summary_combined,
                      list(
                        "winmodel" = winmodel.pkanalysis(obj),
                           "rel_like" = rel_like.pkanalysis(obj)
                           )
                      )

return(summary_combined)
}

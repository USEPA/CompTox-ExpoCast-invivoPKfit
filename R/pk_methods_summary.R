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

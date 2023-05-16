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
#'   all of the models in `obj$settings_optimx$method`.
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
  #ensure that the model has been fitted
  check <- check_required_status(obj = obj,
                                 required_status = status_fit)
  if(!(check %in% TRUE)){
    stop(attr(check, "msg"))
  }

  if(is.null(model)) model <- names(obj$stat_model)
  if(is.null(method)) method <- obj$settings_optimx$method

  sapply(model,
         function(this_model){
           npar <- attr(obj$stat_model[[this_model]]$fit, "npar")
           if(!is.null(npar)){
           fit_par <- obj$stat_model[[this_model]]$fit[method, 1:npar]

           #Add any "constant" params
           if(any(obj$stat_model[[this_model]]$par_DF$optimize_param %in% FALSE &
                  obj$stat_model[[this_model]]$par_DF$use_param %in% TRUE)){
             const_parDF <- subset(obj$stat_model[[this_model]]$par_DF,
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
           }
           }else{
             fit_par <- matrix(data = NA_real_,
                               ncol = sum(obj$stat_model[[this_model]]$par_DF$use_param) +
                                 sum(obj$stat_error_model$sigma_DF$use_param),
                               nrow = length(method))
             rownames(fit_par) <- method
             colnames(fit_par) <- c(
               obj$stat_model[[this_model]]$par_DF[
                 obj$stat_model[[this_model]]$par_DF$use_param %in% TRUE,
                 "param_name"
               ],
               obj$stat_error_model$sigma_DF[
                 obj$stat_error_model$sigma_DF$use_param %in% TRUE,
                 "param_name"
               ]
             )
           }
           return(fit_par)
         },
         USE.NAMES = TRUE,
         simplify = FALSE)
}

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
           }
           return(fit_par)
         },
         USE.NAMES = TRUE,
         simplify = FALSE)
}

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
#'   all of the models in `obj$settings_optimx$method`.
#' @return A named list of numeric matrixes. There is one list element named for
#'   each model in `model`. Each list element is a matrix with as many
#'   rows as items in `method`. The row names are the method names. The matrix column names are
#'   the names of the fitted parameters, including any error standard deviation
#'   hyperparameters (whose names begin with "sigma").
#' @export
#' @author Caroline Ring
#' @family methods for fitted pk objects
coef_sd.pk <- function(obj,
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
           npar <- attr(obj$fit[[this_model]], "npar")
           fit_par <- as.matrix(obj$fit[[this_model]][method, 1:npar])
           #Add any "constant" params
           if(any(obj$prefit[[this_model]]$par_DF$optimize_param %in% FALSE &
                  obj$prefit[[this_model]]$par_DF$use_param %in% TRUE)){
             const_parDF <- subset(obj$prefit[[this_model]]$par_DF,
                                   optimize_param %in% FALSE &
                                     use_param %in% TRUE)[c("param_name",
                                                            "start")]
             const_par <- const_parDF[["start"]]
             names(const_par) <- const_parDF[["param_name"]]
           }else{
             const_par <- NULL
           }

         t(
           apply(fit_par,
                 1,
                 function(this_par){
                   #Evaluate the Hessian
                   numhess <- numDeriv::hessian(func = function(x){
                     log_likelihood(x,
                                    const_params = const_par,
                                    fitdata = obj$data,
                                    data_sigma_group =obj$prefit$stat_error_model$data_sigma_group,
                                    modelfun = obj$stat_model[[this_model]]$conc_fun,
                                    scales_conc = obj$scales$conc,
                                    negative = TRUE,
                                    force_finite = TRUE)
                   },
                   x = this_par,
                   method = 'Richardson')

                   #Invert the Hessian to get the SDs
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
                   names(sds) <- names(this_par)
                   const_sd <- rep(NA_real_, length(const_par))
                   names(const_sd) <- names(const_par)
                   c(sds, const_sd)
                 })
         )
           },
         USE.NAMES = TRUE,
         simplify = FALSE)
}

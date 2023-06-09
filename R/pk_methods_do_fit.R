#' Do fitting
#'
#' Fit PK model(s) for a `pk` object
#'
#' This function estimates the parameters for each model in `stat_model` from
#' the data, using numerical optimization implemented in [optimx::optimx()]. The
#' optimization is done by maximizing the log-likelihood function implemented in
#' [log_likelihood()] (technically, by minimizing the negative log-likelihood).
#' Only the non-excluded observations are used.
#'
#' Due to limitations of [optimx::optimx()], the log-likelihood function is
#' forced to return finite values during this optimization. Impossible
#' combinations of parameters (e.g., parameter values that produce negative
#' predicted concentrations) should have a log-likelihood of `-Inf`, but due to
#' this limitation, they instead have a log-likelihood of `-Machine.doublexmax`.
#' This limitation means that the log-likelihood function is flat in regions of
#' impossible parameter values. It is unlikely, but possible, that the optimizer
#' might get "stuck" in such a flat region -- report convergence, but return a
#' "bad" set of parameter values that produces non-physical predictions.
#'
#' Before trusting the results of any fit, it is recommended to check the
#' log-likelihood using [logLik()] and the Akaike Information Criterion using
#' [AIC()], which check the log-likelihood *without* forcing it to return finite
#' values.
#'
#' @param obj A [pk] object.
#' @return The same [pk] object, with element `fit` containing the fitted
#'   results for each model in `stat_model`.
#' @export
#' @author Caroline Ring
do_fit.pk <- function(obj){
  #check status
  objname <- deparse(substitute(obj))
  status <- obj$status
  if(status >= status_fit){
    warning(paste0(objname,
                   " current status is ",
                   status,
                   ". do_fit() will reset its status to ",
                   status_fit,
                   ". Any results from later workflow stages will be lost."))
  }

  #if preprocessing not already done, do it
  if(obj$status < status_preprocess){
    obj <- do_preprocess(obj)
  }

  #if data summary not already done, do it
  if(obj$status < status_data_info){
    obj <- do_data_info(obj)
  }

  #if prefitting not already done, do it
  if(obj$status < status_prefit){
    obj <- do_prefit(obj)
  }




  #pull the non-excluded observations for fitting
  data <- subset(obj$data, exclude %in% FALSE)
  #pull the sigma parameter corresponding to each non-excluded observation
  data_sigma_group <- data$data_sigma_group

  suppress.messages <- obj$settings_preprocess$suppress.messages
  data_group <- get_data_group(obj)
  data_group_vars <- sapply(data_group,
                            rlang::as_label)

  #For each model:


  fit_list <- sapply(names(obj$stat_model),
         function(this_model){

    if(suppress.messages %in% FALSE){
      message(paste("do_fit.pk(): Fitting model",
                    this_model,
                    "using optimx::optimx()"))
    }

    #nest the necessary data frames...
    data_nest <- get_data(obj) %>%
      tidyr::nest(data = !tidyselect::all_of(data_group_vars))

    par_DF_nest <- obj$prefit$par_DF %>%
      dplyr::filter(model %in% this_model) %>%
      dplyr::select(!model) %>%
      tidyr::nest(par_DF = !tidyselect::all_of(data_group_vars))

    sigma_DF_nest <-  obj$prefit$stat_error_model$sigma_DF %>%
      tidyr::nest(sigma_DF = !tidyselect::all_of(data_group_vars))

    fit_check <- obj$prefit$fit_check %>%
      dplyr::filter(model %in% this_model) %>%
      dplyr::select(!model) %>%
      dplyr::select(!c(n_par, n_sigma, n_detect, n_par_opt, fit_reason))

    #merge it all together
    info_nest <- dplyr::inner_join(
      dplyr::inner_join(
      dplyr::inner_join(data_nest,
                                   par_DF_nest,
                                   by = data_group_vars),
      sigma_DF_nest,
      by = data_group_vars),
      fit_check,
      by = data_group_vars)

    #now use purrr::map2() to run the fit for each group
    fit_out <- info_nest %>%
      dplyr::group_by(!!!data_group) %>%
      dplyr::summarise(fit = {
        if(suppress.messages %in% FALSE){
          cur_data_summary <- dplyr::inner_join(get_data_summary(obj,
                                                                 summary_group = unique(
                                                                   c(data_group,
                                                                     vars(Route,
                                                                          Media)
                                                                   )
                                                                 )),
                                              dplyr::cur_group(),
                                              by = data_group_vars) %>%
          as.data.frame
        message(paste("do_fit.pk(): Fitting model",
                      this_model,
                      "using optimx::optimx()"))
        print(cur_data_summary)
        }

        purrr::pmap(.l = dplyr::cur_data(),
                                   .f = fit_group,
                                   this_model = this_model,
                                   settings_optimx = get_settings_optimx(obj),
                                   modelfun = obj$stat_model[[this_model]]$conc_fun,
                                   dose_norm = obj$scales$conc$dose_norm,
                                   log10_trans = obj$scales$conc$log10_trans,
                                   suppress.messages = suppress.messages)
      }
      )

   fit_out

  }, #end loop over models
  simplify = FALSE,
  USE.NAMES = TRUE)

  obj$fit <- do.call(dplyr::bind_rows,
                     c(fit_list,
                     list(.id = "model")))

  obj$status <- status_fit #fitting complete
  return(obj)
}

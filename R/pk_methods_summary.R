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
  objname <- deparse(substitute(obj))
  status <- get_status(obj)

  cat(paste0(objname, ":", "\n"))
  cat(paste0(attr(status, "msg"), "\n"))

  if(status >= 1){
    #things that require status 1
    #things that are just instructions:
    #data
    #preprocessing settings
    settings_preprocess <- get_settings_preprocess(obj)
    cat(paste0(
      "\nData preprocessing settings:\n",
      paste(names(settings_preprocess),
               settings_preprocess,
               sep = " = ",
               collapse = "\n"),
      "\n"))

    #data info settings
    settings_data_info <- get_settings_data_info(obj)
    cat(paste0(
      "\nData info settings:\n",
      paste(names(settings_data_info),
            settings_data_info,
            sep = " = ",
            collapse = "\n"),
      "\n"))
    #stat_model

models <- names(obj$stat_model)

cat(paste0("\nModels to be fitted:\n",
           paste(models, collapse = ", "),
           "\n"))

    #stat_error_model

error_group <- get_error_group(obj)

cat(paste0("\nError model grouping: ",
    rlang::as_label(error_group),
    "\n"))

    #optimx settings

settings_optimx <- get_settings_optimx(obj)
cat(paste0("\nSettings for optimx::optimx():\n",
           paste(names(settings_optimx),
                 settings_optimx,
                 sep = " = ",
                 collapse = "\n"),
           "\n"))
  }

  #scale_conc

  scale_conc <- get_scale_conc(obj)
  cat(paste0(
    "\nConcentration scaling/transformation:\n",
    paste(names(scale_conc),
          scale_conc,
          sep = " = ",
          collapse = "\n"),
    "\n"))

  #scale_time

  scale_time <- get_scale_time(obj)
  cat(paste0(
    "\nTime scaling/transformation:\n",
    paste(names(scale_time),
          scale_time,
          sep = " = ",
          collapse = "\n"),
    "\n"))

  #Produce a data.frame with all of the instructions in one place
settings_preprocess_DF <- data.frame(instructions_category = "settings_preprocess",
                                     instruction_name = names(settings_preprocess),
                                     instruction_value = sapply(settings_preprocess,
                                                            rlang::as_label))

settings_data_info_DF <- data.frame(instructions_category = "settings_data_info",
                                     instruction_name = names(settings_data_info),
                                     instruction_value = sapply(settings_data_info,
                                                            rlang::as_label))

settings_optimx_DF <- data.frame(instructions_category = "settings_optimx",
                                    instruction_name = names(settings_optimx),
                                    instruction_value = sapply(settings_optimx,
                                                               rlang::as_label))

stat_error_model_DF <- data.frame(instructions_category = "stat_error_model",
                               instruction_name = "error_group",
                               instruction_value = rlang::as_label(error_group))

scale_conc_DF <- data.frame(instructions_category = "scale_conc",
                            instruction_name = names(scale_conc),
                            instruction_value = sapply(scale_conc,
                                                   rlang::as_label))

scale_time_DF <- data.frame(instructions_category = "scale_time",
                            instruction_name = names(scale_time),
                            instruction_value = sapply(scale_time,
                                                   rlang::expr_deparse))

instructions_DF <- do.call(rbind,
                           list(settings_preprocess_DF,
                                settings_data_info_DF,
                                scale_conc_DF,
                                scale_time_DF,
                                stat_error_model_DF,
                                settings_optimx_DF))

rownames(instructions_DF) <- NULL

  #things that require preprocessing & data_info (status 3):
  #data_summary
  cat("Data summary:\n")
  data_info <- get_data_info(obj)
  print(data_info)
  cat("\n")

  #nca
  cat("NCA results:\n")
  nca <- get_nca(obj)
  print(nca)
  cat("\n")

  #things that require pre-fitting (status 4):
#stat_error_model

  #stat_model par_DFs

  #things that require fitting (status 5):
  #coefficients
  #get model coefficients
  coefs <- coef(obj)
  #coef SDs
  #get coefficient SDs
  coef_sds <- coef_sd(obj)
  #goodness of fit metrics

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
                                                      this_outDF <- obj$prefit[[this_model]]$par_DF
                                                      #rowbind sigma_DF
                                                      this_outDF <- rbind(this_outDF,
                                                                          obj$prefit$stat_error_model$sigma_DF)
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

  rownames(outDF) <- NULL

  outDF <- outDF[, c("model", "method", "param_name", "param_units",
                     "optimize_param", "use_param",
                     "lower_bound", "upper_bound", "start",
                     "fitted_value", "fitted_sd")]

  #goodness of fit metrics
  ll_all <- logLik(obj)
  aic_all <- AIC(obj)
  bic_all <- BIC(obj)
  rmse_all <- rmse(obj)
  fold_errors_all <- fold_errors(obj)

  gof_DF_list <- sapply(names(obj$stat_model),
         function(this_model){
           this_ll <- ll_all[[this_model]]
           this_aic <- aic_all[[this_model]]
           this_bic <- bic_all[[this_model]]
           this_rmse <- rmse_all[[this_model]]
           this_fold_errors <- fold_errors_all[[this_model]]
           this_gof_DF <- data.frame(model = this_model,
                                     method = obj$settings_optimx$method)
           this_gof_DF$logLik <- this_ll[this_gof_DF$method]
           this_gof_DF$AIC <- this_aic[this_gof_DF$method]
           this_gof_DF$BIC <- this_bic[this_gof_DF$method]
           this_gof_DF$RMSE <- this_rmse[this_gof_DF$method]
           this_gof_DF$avg_fold_err <- mean(this_fold_errors[, this_gof_DF$method],
                                            na.rm = TRUE)
           this_gof_DF$median_fold_err <- median(this_fold_errors[, this_gof_DF$method],
                                            na.rm = TRUE)
           this_gof_DF$frac_within_2fold <- apply(this_fold_errors[, this_gof_DF$method],
                                                  2,
                                                  function(x) sum(x >= 0.5 & x <= 2)/length(x))
           return(this_gof_DF)
         },
         simplify = FALSE,
         USE.NAMES = FALSE)

  gof_DF <- do.call(rbind, gof_DF_list)

  #model comparison
  model_compare <- compare_models(obj)
  winmodel_DF <- model_compare %>%
    dplyr::group_by(method) %>%
    summarise(winning_model = model[which.min(AIC)]) %>%
    as.data.frame()

  #get TK stats
  tkstats <- get_tkstats(obj)
  #pivot wider and rowbind?


  #compare Cmax, AUC of winning model (for each method) to NCA Cmax, AUC, for each NCA group
  tkwin <- do.call(rbind,
                   mapply(function(this_method, this_model){
    tmp <- get_tkstats(obj, method = this_method, model = this_model)[[1]]
    tmp$model <- this_model
    tmp
  },
  this_method = winmodel_DF$method,
  this_model = winmodel_DF$winning_model,
  SIMPLIFY = FALSE,
  USE.NAMES = FALSE))

  nca <- get_nca(obj)
  tk_nca <- merge(tkwin,
        nca,
        by = setdiff(intersect(names(tkwin),
                               names(nca)),
                     "param_value"),
        suffixes = c(".tkstats", ".nca"))


  return(list("instructions" = instructions_DF,
              "data_summary" = get_data_info(obj),
              "nca" = get_nca(obj),
              "parameters" = outDF,
              "goodness_of_fit" = gof_DF,
              "winning_model" = winmodel_DF,
              "tkstats" = tkstats,
              "tkstats_winning_vs_nca" = tk_nca))
}

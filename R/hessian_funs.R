calc_hessian <- function(pars_opt,
                         pars_const,
                         observations,
                         modelfun,
                         dose_norm,
                         log10_trans){
  hess <- numDeriv::hessian(func = function(x, pars_const) {
    log_likelihood(par = x,
                   const_params = pars_const,
                   data = observations,
                   data_sigma_group = observations$data_sigma_group,
                   modelfun = modelfun,
                   dose_norm = dose_norm,
                   log10_trans = log10_trans,
                   negative = TRUE,
                   force_finite = TRUE)
  },
  x = pars_opt,
  pars_const = pars_const,
  method = 'Richardson')

  colnames(hess) <- names(pars_opt)
  rownames(hess) <- names(pars_opt)

  return(hess)
}


calc_sds_alerts <- function(pars_opt,
                            pars_const,
                            observations,
                            modelfun,
                            dose_norm,
                            log10_trans){

  hess <- calc_hessian(pars_opt = pars_opt,
                       pars_const = pars_const,
                       observations = observations,
                       modelfun = modelfun,
                       dose_norm = dose_norm,
                       log10_trans = log10_trans)

 #  tmp_fun <- function(hess){
 #    diag(solve(hess))^(1 / 2) %>% as.numeric()
 #  }
 #
 #  tmp_fun2 <- function(hess){
 #    V <- Matrix::Cholesky(MASS::ginv(hess),
 #                                         perm = TRUE)
 #
 #                   P1. <- Matrix::expand1(V, which = "P1.")
 #                   P1 <- Matrix::expand1(V, which = "P1")
 #                   L <- Matrix::expand1(V, which = "L")
 #                   L. <- Matrix::expand1(V, which = "L.")
 #
 #                   A <- P1. %*%
 #                     L %*%
 #                     L. %*%
 #                     P1
 #
 #                   diag(as.matrix(A))^(1/2)
 #  }
 #
 #  hess_fun <- safely(
 #      tmp_fun,
 #      otherwise = safely(
 #    quietly(
 #      tmp_fun2
 #    )
 #  )
 #  )
 #
 #  browser()
 # tmp <- safely(tmp_fun)(hess)
 #
 # if(!is.null(tmp$error)){
 #   tmp2 <- safely(
 #     quietly(
 #       tmp_fun2
 #     )
 #   )(hess)
 # }else{
 #   tmp2 <- list(result = NULL,
 #                error = NULL)
 # }
 #
 #
 # output <- if(!is.null(tmp2$))
 #
 #
 #
 #  if(is.null(tmp$error)){
 #    if(is.null(tmp$result$error)){
 #      output <- data.frame(param_sd = set_names(tmp$result$result,
 #                                                names(pars_opt)),
 #                           sd_warning = tmp$result$warnings,
 #                           sd_error = "")
 #    }else{
 #      output <- data.frame(param_sd = set_names(tmp$result$result,
 #                                                names(pars_opt)),
 #                           sd_warning = "",
 #                           sd_error = as.character(tmp$result$error))
 #    }
 #
 #  }else{
 #    output <- data.frame(param_sd = set_names(tmp$result,
 #                                              names(pars_opt)),
 #                         sd_warning = "",
 #                         sd_error = as.character(tmp$error))
 #  }
 #
 #
 #
 #





  output <- tryCatch(
    withCallingHandlers(
    data.frame(param_sd = diag(solve(hess))^(1 / 2) %>% as.numeric(),
         sd_alert = "Hessian inverted using solve()"),
    warning = function(wrn){
      data.frame(param_sd = diag(solve(hess))^(1 / 2) %>% as.numeric(),
                 sd_alert = paste0("Hessian inverted using solve() with warning '",
                                  wrn$message,
                                  "'")
                 )
      invokeRestart("muffleWarning")
    }
    ),
           error = function(err){
             tryCatch(
               withCallingHandlers(
               # pseudovariance matrix
               # see http://gking.harvard.edu/files/help.pdf
               data.frame(param_sd =
                            {
                              V <- Matrix::Cholesky(MASS::ginv(hess),
                                                    perm = TRUE)

                              P1. <- Matrix::expand1(V, which = "P1.")
                              P1 <- Matrix::expand1(V, which = "P1")
                              L <- Matrix::expand1(V, which = "L")
                              L. <- Matrix::expand1(V, which = "L.")

                              A <- P1. %*%
                                L %*%
                                L. %*%
                                P1

                              diag(as.matrix(A))^(1/2)
                            },
                    sd_alert = paste0("solve() failed with error '",
                                      err$message,
                    "'; using pseudovariance matrix ",
                                    "to estimate parameter uncertainty.")
               ),
               warning = function(wrn){
                 data.frame(param_sd =
                              {
                                V <- Matrix::Cholesky(MASS::ginv(hess),
                                                      perm = TRUE)

                                P1. <- Matrix::expand1(V, which = "P1.")
                                P1 <- Matrix::expand1(V, which = "P1")
                                L <- Matrix::expand1(V, which = "L")
                                L. <- Matrix::expand1(V, which = "L.")

                                A <- P1. %*%
                                  L %*%
                                  L. %*%
                                  P1

                                diag(as.matrix(A))^(1/2)
                              },
                            sd_alert = paste0("solve() failed with error '",
                                              err$message,
                                              "'; using pseudovariance matrix ",
                                              "to estimate parameter uncertainty; ",
                                              "warning occurred: '",
                                              wrn$message,
                                              "'")
                 )
                 invokeRestart("muffleWarning")
               }
               ),
               error = function(err2) {
                 data.frame(param_sd = rep(NA_real_, nrow(hess)),
                      sd_alert = paste0("solve() failed with error '",
                      err$message,
                      "' and pseudovariance matrix approach failed with error '",
                      err2$message,
                      "'; returning NAs")
                 )
               }
               )
           }
    )

  output$param_name <- names(pars_opt)

  return(output)
}

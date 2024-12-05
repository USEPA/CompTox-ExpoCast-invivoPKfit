#'Calculate Hessian
#'
#'Calculate Hessian matrix given parameter values and data
#'
#'Calculate the Hessian matrix: the matrix of second derivatives of the
#'objective function with respect to parameters, evaluated for a single set of
#' parameter values for a single model and a single data set.
#' Here, the objective function is the negative
#'log-likelihood implemented in [log_likelihood()], evaluated jointly across the
#'data that was used to fit the model.

#'This is a workhorse function called by [hessian.pk()] and, indirectly, by
#'[coef_sd.pk()].
#'
#'@param pars_opt Named numeric: A vector of parameter values for the parameters
#'  that were optimized. If the length of this vector is \eq{n}, the Hessian
#'  matrix will be \eq{n \times n}. For example, you can get this using
#'  [coef.pk()] with `include_type = "optim"`.
#'@param pars_const Named numeric: A vector of parameter values for parameters
#'  that were held constant, not optimized (but are necessary to evaluate the
#'  model). For example, you can get this using [coef.pk()] with `include_type =
#'  "const"`.
#'@param observations The data used to fit the model. For example, you can get
#'  this using [get_data.pk()].
#'@param modelfun The name of the function that evaluates the model (passed to
#'  [log_likelihood()]).
#'@param dose_norm Logical: Whether to dose-normalize concentrations before
#'  evaluating log-likelihood. Passed to [log_likelihood()].
#'@param log10_trans Logical: Whether to log10-transform concentrations before
#'  evaluating log-likelihood. Passed to [log_likelihood()].
#'
#'@return A square numeric matrix, both dimensions the same as the length of
#'  `pars_opt`. It will have rownames and column names that are the same as the
#'  names of `pars_opt`.
#'
#'@author Caroline Ring
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

#' Inverse diagonal, method 1
#'
#' Get square root of diagonal of inverse matrix, first method
#'
#' Invert a matrix `m` using [solve()], then take the square root of the diagonal.
#'
#' @param m A square numeric matrix, \eq{n \times n}.
#' @return A numeric vector of length \eq{n}.
#'
#' @author Caroline Ring
#'
hess_sd1 <- function(m){
  as.numeric(diag(solve(m)))^(1 / 2)
}

#' Inverse diagonal, method 2
#'
#' Get square root of diagonal of inverse matrix, second method
#'
#' Following the procedure outlined in Gill & King
#' (2004): Calculate generalized inverse of a matrix `m` using [MASS::ginv()]. Then
#' perform a generalized Cholesky factorization of the generalized inverse using
#' [Matrix::Cholesky()] with `perm = TRUE`. Reconstruct the generalized inverse as
#'
#' \deq{ \left( m^{-1} + E \right) = P_1^' L L^' P_1}
#'
#' This should ensure positive semi-definiteness of the reconstruction.
#'
#' Then, take the diagonal of \eq{\left( m^{-1} + E \right)}, and take the square root.
#'
#' @param m A square numeric matrix, \eq{n \times n}.
#' @return A numeric vector of length \eq{n}.
#'
#' @importFrom MASS ginv
#' @importFrom Matrix Cholesky expand1
#'
#' @author Caroline Ring
#'
#' @references Gill J, King G. (2004) What to Do When Your Hessian is Not
#'   Invertible: Alternatives to Model Respecification in Nonlinear Estimation.
#'   Sociological Methods & Research 33(1):54-87. DOI: 10.1177/0049124103262681
#'
hess_sd2 <- function(m){
     m_inv <- MASS::ginv(m)

     #generalized Cholesky factorization
     suppressWarnings(
       V <- Matrix::Cholesky(m_inv, perm = TRUE)
     )

     #reconstruct m_inv from the generalized Cholesky factorization
     #this should ensure positive definiteness
     P1. <- Matrix::expand1(V, which = "P1.")
     P1 <- Matrix::expand1(V, which = "P1")
     L <- Matrix::expand1(V, which = "L")
     L. <- Matrix::expand1(V, which = "L.")

     #A approximates m_inv
     A <- P1. %*%
       L %*%
       L. %*%
       P1

     #have to convert to a standard matrix object
     #in order to take diagonal
     A <- as.matrix(A)

     as.numeric(diag(A))^(1/2)
   }

#' Calculate parameter SDs
#'
#' Calculate parameter SDs using inverse Hessian
#'
#' Calculate parameter SDs using inverse Hessian approach for a single set of
#' parameter values for a single model and a single data set.
#'
#' This is a workhorse function called by [coef_sd.pk()].
#'
#' The coefficient standard deviations are estimated by computing a numerical
#' approximation to the model Hessian (the matrix of second derivatives of the
#' model objective function with respect to each model parameter) and then
#' attempting to invert it. This procedure yields a variance/covariance matrix
#' for the model parameters. The square root of the diagonal elements of this
#' matrix represent the parameter standard deviations.
#'
#' A first attempt is made to invert the Hessian using [solve()] (see
#' [hess_sd1()]). If the Hessian is singular, an attempt is made to calculate a
#' pseudovariance matrix, following the procedure outlined in Gill & King (2004)
#' (see [hess_sd2()]). First, the generalized inverse of the Hessian is
#' calculated using [MASS::ginv()]. Then, a generalized Cholesky decomposition
#' (to ensure positive-definiteness) is calculated using [Matrix::Cholesky] with
#' argument `perm = TRUE`. The generalized inverse is reconstructed from the
#' generalized Cholesky factorization. The square root of the diagonal elements of this matrix
#' represent the parameter standard deviations.
#'
#' If neither of these procedures is successful, then `NA_real_` is returned for
#' all coefficient standard deviations. Record any error messages encountered
#' during the process, and note which method was used to produce the final
#' results. This is a workhorse function called by [coef_sd.pk()].
#'
#' @param pars_opt Named numeric: A vector of parameter values for the
#'   parameters that were optimized. If the length of this vector is \eq{n}, the
#'   Hessian matrix will be \eq{n \times n}. For example, you can get this using
#'   [coef.pk()] with `include_type = "optim"`.
#' @param pars_const Named numeric: A vector of parameter values for parameters
#'   that were held constant, not optimized (but are necessary to evaluate the
#'   model). For example, you can get this using [coef.pk()] with `include_type
#'   = "const"`.
#' @param observations The data used to fit the model. For example, you can get
#'   this using [get_data.pk()].
#' @param modelfun The name of the function that evaluates the model (passed to
#'   [log_likelihood()]).
#' @param dose_norm Logical: Whether to dose-normalize concentrations before
#'   evaluating log-likelihood. Passed to [log_likelihood()].
#' @param log10_trans Logical: Whether to log10-transform concentrations before
#'   evaluating log-likelihood. Passed to [log_likelihood()].
#'
#' @return A data.frame with variables `param_name`, `param_sd`, and `sd_alert`,
#'   and as many rows as the length of `pars_opt`. `param_name` contains the
#'   names of `pars_opt`. `param_sd` contains the parameter standard deviations
#'   calculated using the inverse Hessian. `sd_alerts` is a character variable
#'   noting any errors encountered while attempting to calculate the parameter
#'   SDs.
#' @author Caroline Ring
#' @references Gill J, King G. (2004) What to Do When Your Hessian is Not
#'   Invertible: Alternatives to Model Respecification in Nonlinear Estimation.
#'   Sociological Methods & Research 33(1):54-87. DOI: 10.1177/0049124103262681
#'
calc_sds_alerts <- function(pars_opt,
                            pars_const,
                            observations,
                            modelfun,
                            dose_norm,
                            log10_trans){
#intiialize output
  output <- data.frame(param_name = names(pars_opt),
                       param_sd = rep(NA_real_, length(pars_opt)),
                       sd_alert = rep("", length(pars_opt)))

#calculate Hessian matrix for negative log-likelihood function
  hess <- calc_hessian(pars_opt = pars_opt,
                       pars_const = pars_const,
                       observations = observations,
                       modelfun = modelfun,
                       dose_norm = dose_norm,
                       log10_trans = log10_trans)

  #try inverting Hessian using solve() first
  hess_inv1 <- purrr::safely(hess_sd1)(hess)

  if(is.null(hess_inv1$error)){  #if this succeeded (no error)
    output$param_sd <- hess_inv1$result
    output$sd_alert <- "No error in solve()."
    if(any(!is.finite(hess_inv1$result))){  #if the result had any NaNs
      output$sd_alert <- paste0(output$sd_alert,
                                " solve() method returned NaNs. ",
                                "Trying pseudovariance method.")
      #try pseudovariance method
      hess_inv2 <- purrr::safely(hess_sd2)(hess)
      if(is.null(hess_inv2$error)){ #if pseudovar method succeeded:
        #set output to the result
        output$param_sd <- hess_inv2$result
        output$sd_alert <- paste0(output$sd_alert,
                                  " No error in pseudovariance method. ",
                                  "Returning results of pseudovariance method.")
      }else{ #if pseudovar method errored
        #record error message
        output$sd_alert <- paste0(output$sd_alert,
                                  "Error in pseudovariance method: '",
                                  as.character(hess_inv2$error),
                                  "'. Returning results of solve() method.")
        output$param_sd <- hess_inv1$result
      }
    }else{ #if hess_sd1() succeeded with no NaN
      #then set the output to the result of hess_sd1()
      output$param_sd <- hess_inv1$result
      output$sd_alert <- paste0(output$sd_alert,
                                "Returning results of solve() method.")
    }

  }else{ #if hess_sd1() errored
    #record error message
    output$sd_alert <- paste0("Error in solve(): '",
                             as.character(hess_inv1$error),
                             "'. Trying pseudovariance method.")
    #try pseudovariance method
    hess_inv2 <- purrr::safely(hess_sd2)(hess)
    if(is.null(hess_inv2$error)){ #if inv2 succeeded:
      #set output to the result
      output$param_sd <- hess_inv2$result
      output$sd_alert <- paste0(output$sd_alert,
             " No error in pseudovariance method. ",
             "Returning results of pseudovariance method.")
    }else{ #if pseudovar method errored, too
      #record error message
      output$sd_alert <- paste0(output$sd_alert,
                                " Error in pseudovariance method: '",
                                as.character(hess_inv2$error),
                                "'. Returning NAs.")
      #set output to NAs
      output$param_sd <- rep(NA_real_, nrow(hess))
    }
  }

  return(output)
}

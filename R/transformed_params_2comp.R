#'Transformed parameters for 2-compartment model
#'
#'@param params A named numeric vector of parameters for the 2-compartment
#'  model. Any missing parameters will be filled with `NA_real_`.
#'@param ... Additional arguments (not currently used).
#'@return A named numeric vector of transformed parameters with elements
#'  "alpha", "beta", "A_iv_unit", "B_iv_unit", "A_oral_unit", "B_oral_unit".
#' @author Caroline Ring, Gilberto Padilla Mercado, John Wambaugh
#' @export transformed_params_2comp
#' @family built-in model functions
#' @family 2-compartment model functions
transformed_params_2comp <- function(params,
                                     ...){

  Fgutabs_V1 <- NULL

  list2env(as.list(params), envir = as.environment(-1))

  #see https://www.boomer.org/c/p4/c19/c1902.php
  #for these equations

  alpha_beta_sum <- kelim + k12 + k21
  alpha_beta_prod <- kelim * k21

  alpha <- (alpha_beta_sum + sqrt(alpha_beta_sum^2 - 4*alpha_beta_prod)) / 2
  beta <- (alpha_beta_sum - sqrt(alpha_beta_sum^2 - 4*alpha_beta_prod)) / 2

  A_iv_unit <- (alpha - k21) /
    (V1 * (alpha - beta))
  B_iv_unit <- (k21 - beta) /
    (V1 * (alpha - beta))

  A_oral_unit <- (kgutabs * Fgutabs_V1 *
                    (alpha - k21)) /
    ( (kgutabs - alpha) * (alpha - beta))

  B_oral_unit <- (kgutabs * Fgutabs_V1 *
                    (k21 - beta)) /
    ( (kgutabs - beta) * (alpha - beta))

  return(c("alpha" = alpha,
              "beta" = beta,
              "A_iv_unit" = A_iv_unit,
              "B_iv_unit" = B_iv_unit,
              "A_oral_unit" = A_oral_unit,
              "B_oral_unit" = B_oral_unit))

}

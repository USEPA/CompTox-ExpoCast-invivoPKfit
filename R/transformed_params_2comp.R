transformed_params_2comp <- function(params){

  params <- fill_params_2comp(params)

  #see https://www.boomer.org/c/p4/c19/c1902.php
  #for these equations

  alpha_beta_sum <- params$kelim + params$k12 + params$k21
  alpha_beta_prod <- params$kelim * params$k21

  alpha <- (alpha_beta_sum + sqrt(alpha_beta_sum^2 - 4*alpha_beta_prod)) / 2
  beta <- (alpha_beta_sum - sqrt(alpha_beta_sum^2 - 4*alpha_beta_prod)) / 2

  A_iv_unit <- (alpha - params$k21) /
    (params$V1 * (alpha - beta))
  B_iv_unit <- (params$k21 - beta) /
    (params$V1 * (alpha - beta))

  A_oral_unit <- (params$kgutabs * params$Fgutabs_V1 *
                    (alpha - params$k21)) /
    ( (params$kgutabs - alpha) * (alpha - beta))

  B_oral_unit <- (params$kgutabs * params$Fgutabs_V1 *
                    (params$k21 - beta)) /
    ( (params$kgutabs - beta) * (alpha - beta))

  return(list("alpha" = alpha,
              "beta" = beta,
              "A_iv_unit" = A_iv_unit,
              "B_iv_unit" = B_iv_unit,
              "A_oral_unit" = A_oral_unit,
              "B_oral_unit" = B_oral_unit))

}

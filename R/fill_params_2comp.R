#' Fill parameters for 2-compartment model
#'
#' @param params Named numeric vector of parameters for the 2-compartment model
#' @return A named numeric vector of parameters, with any 2-compartment model
#'   parameters not present in `params` filled with `NA_real_`. If any two of
#'   `Fgutabs`, `V1`, and `Fgutabs_V1` were present in `params`, the third will
#'   be imputed to agree with the other two.
#' @author Caroline Ring

fill_params_2comp <- function(params){

  #fill in missing params with NAs
  missing_params <- setdiff(model_2comp$params,
                            names(params))
  params[missing_params] <- NA_real_

  if(is.na(params["Fgutabs_V1"])){
    params["Fgutabs_V1"] <- params["Fgutabs"]/params["V1"]
  }

 if(is.na(params["Fgutabs"])){
    params["Fgutabs"] <- params["Fgutabs_V1"] * params["V1"]
 }

if(is.na(params["V1"])){
    params["V1"] <- params["Fgutabs"] / params["Fgutabs_V1"]
}

  return(params)
}

#' Fill parameters for 1-compartment model
#'
#' @param params Named numeric vector of parameters for the 1-compartment model
#' @return A named numeric vector of parameters, with any 1-compartment model
#'   parameters not present in `params` filled with `NA_real_`. If any two of
#'   `Fgutabs`, `Vdist`, and `Fgutabs_Vdist` were present in `params`, the third will
#'   be imputed to agree with the other two.
#' @author Caroline Ring

fill_params_1comp <- function(params) {

  # fill in missing params with NAs
  missing_params <- base::setdiff(model_1comp$params,
                            names(params))
  params[missing_params] <- NA_real_

  if (is.na(params["Fgutabs_Vdist"])) {
    params["Fgutabs_Vdist"] <- params["Fgutabs"] / params["Vdist"]
  }

 if (is.na(params["Fgutabs"])) {
    params["Fgutabs"] <- params["Fgutabs_Vdist"] * params["Vdist"]
 }

if (is.na(params["Vdist"])) {
    params["Vdist"] <- params["Fgutabs"] / params["Fgutabs_Vdist"]
}

  return(params)
}

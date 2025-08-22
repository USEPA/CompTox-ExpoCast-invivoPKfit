#' Fill parameters for 2-compartment model
#'
#' @param params Named list of parameters for the 2-compartment model.
#' @return A named numeric list of parameters, with any 2-compartment model
#'   parameters not present in `params` filled with `NA_real_`. If any two of
#'   `Fgutabs`, `V1`, and `Fgutabs_V1` were present in `params`, the third will
#'   be imputed to agree with the other two.
#' @author Caroline Ring

fill_params_2comp <- function(params) {

  # fill in missing params with NAs
  missing_params <- base::setdiff(model_2comp$params, names(params))
  params[missing_params] <- NA_real_

  has_Fgutabs <- !"Fgutabs" %in% missing_params
  has_V1 <- !"V1" %in% missing_params
  has_Fgutabs_V1 <- !"Fgutabs_V1" %in% missing_params

  if (!has_Fgutabs_V1 && has_Fgutabs && has_V1) {
    params["Fgutabs_V1"] <- params[["Fgutabs"]] / params[["V1"]]
  }

  if (!has_Fgutabs && has_Fgutabs_V1 && has_V1) {
    params["Fgutabs"] <- params[["Fgutabs_V1"]] * params[["V1"]]
  }

  if (!has_V1 && has_Fgutabs_V1 && has_Fgutabs) {
    params["V1"] <- params[["Fgutabs"]] / params[["Fgutabs_V1"]]
  }

  return(params)
}

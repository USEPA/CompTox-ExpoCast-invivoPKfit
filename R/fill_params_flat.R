#' Fill parameters for flat model
#'
#' @param params Named list of parameters for the flat model.
#' @return A named numeric vector of parameters, with any flat model
#'   parameters not present in `params` filled with `NA_real_`. If any two of
#'   `Fgutabs`, `Vdist`, and `Fgutabs_Vdist` were present in `params`, the third will
#'   be imputed to agree with the other two.
#' @author Caroline Ring

fill_params_flat <- function(params) {

  # fill in missing params with NAs
  missing_params <- base::setdiff(model_flat$params,
                            names(params))
  params[missing_params] <- NA_real_

  has_Fgutabs <- !"Fgutabs" %in% missing_params
  has_Vdist <- !"Vdist" %in% missing_params
  has_Fgutabs_Vdist <- !"Fgutabs_Vdist" %in% missing_params

  if (!has_Fgutabs_Vdist && has_Fgutabs && has_Vdist) {
    params["Fgutabs_Vdist"] <- params[["Fgutabs"]] / params[["Vdist"]]
  }

  if (!has_Fgutabs && has_Fgutabs_Vdist && has_Vdist) {
    params["Fgutabs"] <- params[["Fgutabs_Vdist"]] * params[["Vdist"]]
  }

  if (!has_Vdist && has_Fgutabs_Vdist && has_Fgutabs) {
    params["Vdist"] <- params[["Fgutabs"]] / params[["Fgutabs_Vdist"]]
  }

  return(params)
}

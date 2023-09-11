#' Check flat model parameters
#'
#' Check to make sure required parameters are present to evaluate flat
#' model for a given route and medium
#'
#' @param params A named numeric vector of parameters for the flat model.
#' @param route A character vector of routes: "iv" and/or "oral".
#' @param medium A character vector of tissue media: "plasma" and/or "blood".
#' @param ... Additional arguments (not currently used)
#' @return Character: A message. If all required parameters are present for the
#'   given media & routes, the message is "Parameters OK". If required
#'   parameters for the oral route are missing, the message is "Error: For
#'   flat oral model, missing parameters (comma-separated list of
#'   parameter names)".  If required
#'   parameters for the IV route are missing, the message is "Error: For
#'   flat oral model, missing parameters (comma-separated list of
#'   parameter names)".
#' @author Caroline Ring

check_params_flat <- function(params,
                               route,
                               medium,
                               ...){

  msg <- "Parameters OK"

  # params <- fill_params_flat(params)

  #check for any missing parameters
  #required params for oral dose
  if(any(route %in% "oral")){
    missing_params <- setdiff(c("Fgutabs_Vdist"),
                              names(params)[is.finite(params)])
    if(length(missing_params)>0){
      msg <- (paste("Error: For flat oral model,",
                 "missing parameters:",
                 paste(missing_params, collapse = ", ")))
    }
  }

  #required params for IV dose
  if(any(route %in% "iv")){
    missing_params <- setdiff(c("Vdist"),
                              names(params)[is.finite(params)])
    if(length(missing_params)>0){
      msg <- (paste("Error: For flat IV model,",
                 "missing parameters:",
                 paste(missing_params, collapse = ", ")))
    }
  }

  if(any(medium %in% "blood")){
    if(!("Rblood2plasma" %in% names(params)[is.finite(params)])){
      msg <- (paste0("Error: For flat model ",
                  "in blood: missing parameter Rblood2plasma"))
    }
  }

  return(msg)
}

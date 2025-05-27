#' Check 1-compartment model parameters with specific clearance
#'
#' Check to make sure required parameters are present to evaluate 1-compartment
#' model for a given route and medium
#'
#' @param params A named numeric vector of parameters for the 1-compartment model.
#' @param route A character vector of routes: "iv" and/or "oral".
#' @param medium A character vector of tissue media: "plasma" and/or "blood".
#' @param ... Additional arguments (not currently used)
#' @return Character: A message. If all required parameters are present for the
#'   given media & routes, the message is "Parameters OK". If required
#'   parameters for the oral route are missing, the message is "Error: For
#'   1-compartment oral model, missing parameters (comma-separated list of
#'   parameter names)".  If required
#'   parameters for the IV route are missing, the message is "Error: For
#'   1-compartment oral model, missing parameters (comma-separated list of
#'   parameter names)".
#' @author Caroline Ring

check_params_1comp_cl <- function(params,
                               route,
                               medium,
                               ...) {

  msg <- "Parameters OK"

  # check for any missing parameters
  # required params for oral dose
  if (any(route %in% "oral")) {
    missing_params <- setdiff(c("Clint",
                                "Q_gfr",
                                "Q_totli",
                                "Q_alv",
                                "Kblood2air",
                                "Fup",
                                "Fgutabs_Vdist",
                                "Rblood2plasma",
                                "kgutabs"),
                              names(params)[is.finite(params)])
    if (length(missing_params) > 0) {
      msg <- (paste("Error: For 1-compartment oral model,",
                    "missing parameters:",
                    toString(missing_params)
                    )
              )
    }
  }

  # required params for IV dose
  if (any(route %in% "iv")) {
    missing_params <- setdiff(c("Clint",
                                "Q_gfr",
                                "Q_totli",
                                "Q_alv",
                                "Kblood2air",
                                "Rblood2plasma",
                                "Fup",
                                "Vdist"),
                              names(params)[is.finite(params)])
    if (length(missing_params) > 0) {
      msg <- (paste("Error: For 1-compartment IV model,",
                    "missing parameters:",
                    toString(missing_params)
                    )
              )
    }
  }

  return(msg)
}

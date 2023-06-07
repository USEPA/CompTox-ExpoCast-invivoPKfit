fill_params_2comp <- function(params){
  if(all(c("Fgutabs", "V1") %in% names(params))){
    params$Fgutabs_V1 <- params$Fgutabs/params$V1
  }

  #if V1 and Fgutabs_V1 are provided, but not Fgutabs, compute Fgutabs
  if(all(c("V1", "Fgutabs_V1") %in% names(params)) &
     !("Fgutabs" %in% names(params))
  ){
    params$Fgutabs <- params$Fgutabs_V1 * params$V1
  }

  #if Fgutabs and Fgutabs_V1 provided, but not V1, compute V1
  if(all(c("Fgutabs", "Fgutabs_V1") %in% names(params)) &
     !("V1" %in% names(params))
  ){
    params$V1 <- params$Fgutabs / params$Fgutabs_V1
  }

  #drop any length-0 params
  param_length <- sapply(params, length)
  params <- params[param_length>0]

  return(params)
}

check_params_2comp <- function(params,
                               route,
                               medium,
                               ...){

  msg <- "Parameters OK"

  params <- fix_params_2comp(params)

  #check for any missing parameters
  #required params for oral dose
  if(any(route %in% "oral")){
    missing_params <- setdiff(c("kelim",
                                "k21",
                                "k12",
                                "Fgutabs_V1",
                                "kgutabs"),
                              names(params))
    if(length(missing_params)>0){
      msg <- (paste("Error: For 2-compartment oral model,",
                 "missing parameters:",
                 paste(missing_params, collapse = ", ")))
    }
  }

  #required params for IV dose
  if(any(route %in% "iv")){
    missing_params <- setdiff(c("kelim",
                                "V1",
                                "k12",
                                "k21"),
                              names(params))
    if(length(missing_params)>0){
      msg <- (paste("Error: For 2-compartment IV model,",
                 "missing parameters:",
                 paste(missing_params, collapse = ", ")))
    }
  }

  if(any(medium %in% "blood")){
    if(!("Rblood2plasma" %in% names(params))){
      msg <- (paste0("Error: For 2-compartment model ",
                  "in blood: missing parameter Rblood2plasma"))
    }
  }

  return(msg)
}

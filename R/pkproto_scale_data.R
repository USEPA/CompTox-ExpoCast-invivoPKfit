#'Scale concentrations
#'
#'Normalize and/or transform concentration data
#'
#'# Requirements for `expr` `expr` should be an expression (unquoted) --
#'effectively, the "right hand side" of a formula that will produce the
#'transformed concentration variable. This expression may refer to any of the
#'harmonized `invivopkfit` variable names (see [pk()] for a list), except that
#'the concentration variable to be transformed (e.g. `Conc`, `Value`,
#'`Value_SD`, etc.) should be replaced with the placeholder name `.conc`. For
#'example, to apply dose-normalization (i.e., new concentration = old
#'concentration / dose), you would supply
#'
#'`expr = .conc/Dose`
#'
#'To apply a log10 transformation (i.e., new concentration = log10(old
#'concentration)), you would supply
#'
#'`expr = log10(.conc)`
#'
#'
#'The reason for using the `.conc` placeholder is that this formula will be
#'applied to five different concentration variables in the harmonized data
#'(`Value`, `Value_SD`, `LOQ`, `Conc`, and `Conc_SD`), and (during curve
#'fitting) will also be applied to model-predicted concentrations. The `.conc`
#'placeholder is an easy way to say "Apply this transformation to all of the
#'concentration-related variables."
#'
#'The expression may additionally refer to variables defined in the global
#'environment, or the environment where you are calling this function. Variables
#'will be looked for in the harmonized `invivopkfit` data frame first, and if
#'not found there, will be looked for in the environment. This means if you have
#'a variable in your global environment that is named the same thing as one of
#'the `invivopkfit` harmonized variables, you should rename your global-environment variable before
#'trying to refer to it in `expr`.
#'
#'For example, if you write at the R command line
#' ```
#' b <- 5
#' my_pk <- pk(data = my_df) + scale_conc(expr = .conc/b) + stat_model(model = "1comp")
#' my_pk <- fit(my_pk)
#'```
#'Then all concentrations will be divided by 5, which you will be able to see by
#'examining `my_pk$data`.
#'
#'
#'@param expr An expression to compute normalized and/or transformed
#'  concentration data. See Details for requirements.
#'@param ... Other arguments (not currently used)
#'@return An object of class `pk_scales`: A named list with element `name =
#'  "conc"` (denoting the variable to be scaled) and `value = list("expr" =
#'  rlang::enquo(expr), ...)` (denoting the arguments supplied to
#'  [scale_conc()]). See [pk_add.pk_scales()].
#'@export
#'@author Caroline Ring
scale_conc <- function(dose_norm = FALSE,
                       log10_trans = FALSE,
                       ratio_conc_dose = 1,
                       ...){

  #Initialize scale_conc as a list with element "name"
  scale_conc <- list(name = "conc",
                     value = list(
                      "ratio_conc_dose" = ratio_conc_dose,
                     "dose_norm" = dose_norm,
                     "log10_trans" = log10_trans)
  )

  #Set up expression
  expr <- quote(.conc)
  #add ratio_conc_dose normalization step to expression
  if(!(ratio_conc_dose %in% 1)){
  expr <- substitute(x/ratio_conc_dose,
                     list(x = expr,
                          ratio_conc_dose = ratio_conc_dose))
  }

  if(dose_norm %in% TRUE){
    #add dose normalization step to expression
    expr <- substitute(x/Dose,
                       list(x = expr))
  }


  if(log10_trans %in% TRUE){
    #add log transformation step to expression
    expr <- substitute(log10(x),
                       list(x = expr))
  }

  #enquote expr
  scale_conc$value$expr <- rlang::new_quosure(expr = expr,
                                              env = rlang::caller_env())
  #this qusoure will have environment "global"

  #get any other arguments and values
  scale_conc$value <- c(scale_conc$value, list(...))
  #set class
  class(scale_conc) <- c(class(scale_conc), "pkproto", "pk_scales")

return(scale_conc)
}

#' Scale times
#'
#' Transform time data
#'
#' @param new_units New units to use for time. Default is `"identity"` (leave
#'   time in the original units). Another useful option is `"auto"`, to
#'   automatically select new time units based on the time of the last detected
#'   observation. You may also specify any time units understood by
#'   `lubridate::duration()`, i.e., `"seconds"`, `"hours"`, `"days"`, `"weeks"`,
#'   `"months"`, `"years"`, `"milliseconds"`, `"microseconds"`, `"nanoseconds"`,
#'   and/or `"picoseconds"`. You may only specify one new unit (e.g., `new_units
#'   = c("days", "weeks")` is not valid).
#' @param ... Other arguments (not currently used)
#' @return An object of class `pk_scales`: A named list with two elements `name
#'   = "time"` (denoting the variable to be scaled) and `value = list("new_units" =
#'   new_units, ...)` (denoting the arguments supplied to [scale_time()]). See
#'   [pk_add.pk_scales()].
scale_time <- function(new_units = "identity",
                          ...){
  #get arguments and values
  argg <- c(as.list(environment()), list(...))
  scale_time <- list(name = "time")
  scale_time$value <- argg
  #set class
  class(scale_time) <- c(class(scale_time), "pkproto", "pk_scales")

  return(scale_time)
}





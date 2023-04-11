#' Scale concentrations
#'
#' Normalize and/or transform concentration data
#'
#' @param normalize A normalization to apply to concentrations. Built-in options
#'   are `"identity"` (the default; i.e., divide by 1) and `"dose_normalize"`
#'   (i.e., divide each concentration value by its corresponding dose; see
#'   [dose_normalize()]). May also be the name of a function. If no function is
#'   found by the provided name in the calling environment, then
#'   [scale_conc.pk()] will attempt to interpret it as a transformation from the
#'   `scales` package. For example, if you provide `"boxcox"`, it will be
#'   interpreted as [scales::boxcox_trans()].
#' @param trans A transformation to apply to concentrations. Will be applied
#'   *after* the normalization. Default is `"identity"` (no transformation). You
#'   can supply the name of your own function. If no function by that name is
#'   found, then this function will try to interpret it as one of the
#'   transformations in the `scales` package.
#' @return An object of class `pk_scales`: A named list with element `variable =
#'   "conc"` (denoting the variable to be scaled) and `options = list(normalize
#'   = normalize, trans = trans)` (denoting the `normalize` and `trans`
#'   arguments supplied to [scale_conc()]). See [pk_add.pk_scales()].
#' @export
#' @author Caroline Ring
scale_conc <- function(normalize = "identity",
                       trans = "identity"){
  scale_conc <- list(variable = "conc")
  scale_conc$options <- list("normalize" = normalize,
                             "trans" = trans)
  class(scale_conc) <- c(class(scale_conc), "pkproto", "pk_scales")

return(scale_conc)
}

#' Scale times
#'
#' Transform time data
#'
#' @param trans A transformation to apply to time data -- *i.e.*, new units to
#'   use for time. Default is `"identity"` (no transformation, leave time in the
#'   original units). Another useful option is `"auto"`, to automatically select
#'   new time units based on the time of the last detected observation. You may
#'   also specify any time units understood by `lubridate::duration()`, i.e.
#'   "seconds", "hours", "days", "weeks", "months", "years", "milliseconds",
#'   "microseconds", "nanoseconds", and/or "picoseconds".
#' @return An object of class `pk_scales`: A named list with two elements
#'   `variable = "time"` (denoting the variable to be scaled) and `options =
#'   list(trans = trans)` (denoting the `normalize` and `trans` arguments
#'   supplied to [scale_time()]). See [pk_add.pk_scales()].
scale_time.pk <- function(obj,
                          trans = "identity"){
  scale_time <- list(variable = "time")
  scale_time$options <- list("trans" = trans)
  class(scale_time) <- c(class(scale_time), "pkproto", "pk_scales")

  return(scale_time)
}

#' Dose-normalize concentration data
#'
#' Apply dose normalization to concentration data.
#'
#' Dose normalization simply divides a tissue concentration (mg/L) by the
#' corresponding administered dose (mg/kg).
#'
#' @param conc A numeric vector of concentrations
#' @param dose The corresponding numeric vector of doses
#' @return A numeric vector, `conc/dose`
#' @author Caroline Ring
dose_normalize <- function(conc, dose){
  conc/dose
}



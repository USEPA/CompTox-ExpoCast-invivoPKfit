#' Check required status
#'
#' Check whether a [pk()] object has a particular required status level
#'
#' This is a helper function to check whether a [pk()] object has the status
#' required for certain operations. For example, status 4 (fitting complete) is
#' required for any fit evaluation functions: [predict.pk()], [residuals.pk()],
#' [coef.pk()], [coef_sd.pk()], [rmse.pk()], [fold_errors.pk()]
#'
#' @param obj A [pk()] object
#' @param required_status Integer: The required status. 1 = initialized; 2 =
#'   pre-processed; 3 = pre-fitted; 4 = fitted.
#' @return If the [pk()] object has the required status or greater, returns
#'   TRUE. If the [pk()] object has less than the required status, returns
#'   FALSE. Returned value has an attribute `msg`, containing an informative
#'   message as a string.
#' @export
#' @author Caroline Ring
check_required_status.pk <- function(obj,
                                   required_status){
  objname <- deparse(substitute(obj))
  status <- obj$status
  if(status > 4) status <- 5
  steps <- c("1/4. Object has been initialized (pk())",
             "2/4. Data pre-processing complete (preprocess_data())",
             "3/4. Model pre-fitting complete (prefit())",
             "4/4. Model fitting complete (fit())",
             "status > 4 (not implemented)")
if(status < required_status){
msg <- paste(
    paste(objname, "current status is below required status."),
paste0("Current status for ", objname, ":"),
steps[status],
"Required status:",
steps[required_status],
sep = "\n")

out <- FALSE
attr(out, "msg") <- msg
}else{
  out <- TRUE
  attr(out, "msg") <- NULL
}
  return(out)
}

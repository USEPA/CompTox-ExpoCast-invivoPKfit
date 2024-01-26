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
#'
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
  if(status > status_fit) status <- status_fit + 1
  steps <- c("1/5. Object has been initialized",
             "2/5. Data pre-processing complete",
             "3/5. Data information summary and NCA complete",
             "4/5. Model pre-fitting complete",
             "5/5. Model fitting complete",
             "status > 5. not implemented")
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

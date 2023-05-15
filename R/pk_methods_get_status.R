#' Check status of a `pk` object
#'
#' Check status of a `pk` object
#'
#' `pk` objects have integer statuses reflecting what stage of the analysis
#' process they are at.
#'
#' 1. Object has been initialized
#' 2. Data pre-processing complete
#' 3. Model pre-fitting complete 4
#' . Model fitting complete
#'
#' If a `pk` object of status 2 or greater has its instructions modified with
#' `+`, then its status will be reset to 1, indicating that any analysis results
#' contained in the object are now outdated and all steps of the analysis need
#' to be re-run.
#'
#' This function allows the user to check the status of a `pk` object.
#'
#' A message will be printed listing the analysis steps that have been completed
#' for this `pk` object, and the integer status will be returned.
#'
#' @param obj A `pk` object
#' @return The status of the `pk` object as an integer.
#' @export
#' @author Caroline Ring
get_status.pk <- function(obj){
  objname <- deparse(substitute(obj))
  obj_status <- obj$status
  steps <- c("1/5. Object has been initialized",
             "2/5. Data pre-processing",
             "3/5. Data information summary and NCA",
             "4/5. Model pre-fitting ",
             "5/5. Model fitting")
  completed_steps <- copy(steps)
  completed_steps[seq(1, obj_status)] <- paste(completed_steps[seq(1, obj_status)],
                                               "--",
                                               "COMPLETE")

  completed_steps[seq(obj_status+1, length(steps))] <- paste(completed_steps[seq(obj_status+1, length(steps))],
                                                           "--",
                                                           "NOT COMPLETE")
  msg <- paste(
    paste0("Status of pk object ", objname, ":"),
    paste(completed_steps,
          collapse = "\n"),
    sep  = "\n"
  )
  message(msg)
  out <- obj_status
  attr(out, "msg") <- msg
  return(obj_status)
}

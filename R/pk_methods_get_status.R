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
get_status.pk <- function(obj){
  steps <- c("1. Object has been initialized",
             "2. Data pre-processing complete",
             "3. Model pre-fitting complete",
             "4. Model fitting complete")
  message(paste(steps[seq(1, obj$status)],
                collapse = "/n"))
  return(obj$status)
}

#' Check new data
#'
#' Check new data to ensure it has the required variables and classes
#'
#' This is a helper function to check new data to ensure it has the required
#' variables and that those variables are of the correct classes. This is
#' useful, for example, when making predictions from a fitted [pk()] model
#' object on new data.
#'
#' @param newdata A `data.frame` containing new data
#' @param olddata A `data.frame` containing existing data. `newdata` variable
#'   classes will be required to match `olddata`
#' @param req_vars A `character` vector of required variable names that must
#'   appear in `newdata`
#' @param exclude Logical: Whether a variable `"exclude"` also must be present in `newdata`
#' @return `TRUE`, if required variables are present in `newdata`, and required
#'   variables are of the same class in `newdata` and `olddata`. Otherwise, this
#'   function will stop with an error.
#' @author Caroline Ring

check_newdata <- function(newdata,
                          olddata,
                          req_vars,
                          exclude = FALSE){

  if(exclude %in% TRUE){
    req_vars <- union(req_vars,
                      "exclude")
  }

  #check that newdata has the required variables
  if(!(all(req_vars %in% names(newdata)))){
    stop("newdata is missing one or more required variables.\n",
         "Required variables: ",
         toString(req_vars), "\n",
         "Missing required variables: ",
         toString(setdiff(req_vars, names(newdata)))
    )
  }

 #check that required variables are the same type as old data

  req_class <- sapply(req_vars,
                     function(this_var) {
                       class(olddata[[this_var]])
                     },
                     simplify = FALSE,
                     USE.NAMES = TRUE
                       )

  new_class <- sapply(req_vars,
                     function(this_var) {
                       class(newdata[[this_var]])
                     },
                     simplify = FALSE,
                     USE.NAMES = TRUE
  )

  has_class <- sapply(seq_along(req_class),
                      function(i){
                        req_class[[i]] %in% new_class[[i]]
                      },
                      simplify = TRUE)

  if(!all(has_class %in% TRUE)){
    bad_class <- req_vars[has_class %in% FALSE]

    stop(paste("The following variables in newdata have the wrong class:",
               paste(bad_class,
                     "-- Required class:",
                     req_class[has_class %in% FALSE],
                     "; newdata class:",
                     new_class[has_class %in% FALSE],
                     sep = " ",
                     collapse = "\n"),
               sep = "\n"
               ))
  }
 #if everything is OK, return TRUE.
  return(TRUE)
}

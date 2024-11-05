#' Subtract a pkproto object from a pk object
#'
#' @param e1 A pk pbject
#' @param e2 A pkproto object
#' @return The pk object, modified by adding the pkproto object
#'@export
#' @author Caroline Ring
"-.pk" <- function(e1, e2) {
  if (missing(e2)) {
    cli::cli_abort(c(
      "Cannot use {.code +} with a single argument",
      "i" = "Did you accidentally put {.code -} on a new line?"
    ))
  }

  # Get the name of what was passed in as e2, and pass along so that it
  # can be displayed in error messages
  e2name <- deparse(substitute(e2))

  if      (is.pk(e1))  subtract_pk(e1, e2, e2name)
  else if (is.pkproto(e1)) {
    cli::cli_abort(c(
      "Cannot subtract {.cls pkproto} objects from each other",
      "i" = "Did you forget to subtract this object from a {.cls pk} object?"
    ))
  }
}

#' Subtract various `pkproto` objects from a `pk` object
#'
#'@param pk_obj The `pk` object
#'@param object The `pkproto` object to be subtracted
#'@param objectname The name of the `pkproto` object to be subtracted
#'
#'@return The `pk` object modified by the subtraction.
#' @import cli
#' @export
subtract_pk <- function(pk_obj, object, objectname) {
  if (is.null(object)) return(pk_obj)

  p <- pk_subtract(object, pk_obj, objectname)
  p
}

#' Subtract a `pk_stat_model` object.
#'
#' @param pkproto_obj The `pk_stat_model` object to be subtracted.
#' @param pk_obj The `pk` object from which the `pk_stat_model` object will be subtracted.
#' @param objectname The name of the `pk_stat_model` object.
#'
#' @return The `pk` object, modified by subtracting the `stat_model`.
#' @author Caroline Ring
#' @export
pk_subtract.pk_stat_model <- function(pkproto_obj, pk_obj, objectname){

  for (this_model in names(pkproto_obj)) {
    if (!is.null(pk_obj$stat_model[[this_model]])) {
      #if this_model is present, remove it, and say so
      message(objectname,
              ": stat_model for",
              this_model,
              "will be removed"
      )
      pk_obj$stat_model[[this_model]] <- NULL

    } else {
      #if this_model not present, say so
      message(objectname,
              ": stat_model for",
              this_model,
              "was not present to begin with, so it can't be removed"
      )
    }
  }
  if (pk_obj$status > 2L) {
    message(objectname,
            ": Modifying stat_model resets status to level 2 ",
            "(data preprocessing complete);\n",
            "model pre-fit (level 3)  (prefit()) ",
            "and model fit (level 4) (fit()) will need to be re-done"
    )
    pk_obj$status <- 2L #data pre-processing won't change with addition of new model, but model pre-fit and fit will change
  }

  return(pk_obj)
}

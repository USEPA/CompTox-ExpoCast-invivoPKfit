#' Add a `pkproto` object to a `pk_faceted` object
#'
#' Add a `pkproto` object to a `pk_faceted` object
#'
#' Note that `e1 + e2` is equivalent to
#'
#' ```
#' `+`(e1, e2)
#' ```
#'
#' @param e1 A `pk_faceted` pbject
#' @param e2 A `pkproto` object
#' @return The `pk` object, modified by adding the `pkproto` object
#'@export
#' @author Caroline Ring
"+.pk_faceted" <- function(e1, e2) {
  #throw error if only one argument
  if (missing(e2)) {
    cli::cli_abort(c(
      "Cannot use {.code +} with a single argument",
      "i" = "Did you accidentally put {.code +} on a new line?"
    ))
  }

  # Get the name of what was passed in as e2, and pass along so that it
  # can be displayed in error messages
  e2name <- deparse(substitute(e2))

  #make sue
  if      (is.pk_faceted(e1)) {
    if(!inherits(e2, "pk_facet_data")){
    #add to each of the pk objects in turn
    e1 %>% dplyr::transmute(pk_object = purrr::map(pk_object,
                                                   function(x){
                                                     add_pk(x, e2, e2name)
                                                   }))
    }else{ #can't re-facet already-faceted data.
      e1name <- deparse(substitute(e1))
      stop(paste(e1name,
                 "is already faceted. New faceting",
                 e2name,
                 "cannot be applied."))
    }
  }
  else if (is.pkproto(e1)) {
    cli::cli_abort(c(
      "Cannot add {.cls pkproto} objects together",
      "i" = "Did you forget to add this object to a {.cls pk} object?"
    ))
  }
}

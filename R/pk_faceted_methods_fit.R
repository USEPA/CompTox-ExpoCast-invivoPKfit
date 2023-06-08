#' Fitting
#'
#' Fit PK model(s) for a `pk_faceted` object
#'
#' @param obj A `pk_faceted` object.
#' @return The same `pk_faceted` object, with list column `pk_object` modified
#'   by applying [fit.pk()] to each of its elements.
#' @export
#' @author Caroline Ring
fit.pk_faceted <- function(obj){
  obj <- obj %>%
    dplyr::transmute(pk_object = purrr::map(pk_object,
                                            fit)
    )

  return(obj)
}

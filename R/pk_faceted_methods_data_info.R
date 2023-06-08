#' calculate summary data info for a `pk_faceted` object
#'
#' Calculate summary data information, including non-compartmental analysis.
#'
#'
#' @param obj An object of class `pk_faceted``
#' @export
#' @importFrom magrittr `%>%`
#' @author Caroline Ring
data_info.pk_faceted <- function(obj){

  facets <- attr(obj, "facets")

  obj <- obj %>%
    dplyr::transmute(pk_object = purrr::map(pk_object,
                                            data_info)
    )
  return(obj)

  return(obj)
}

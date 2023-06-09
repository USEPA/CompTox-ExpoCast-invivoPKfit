#' Pre-process data
#'
#' Pre-process data for a `pk_faceted` object
#'
#' For each group of data defined by unique combinations of the faceting
#' variables,
#'
#'
#' @param obj A `pk_faceted` object
#' @return The same `pk_faceted` object, with list column `pk_object` modified
#'   by applying [preprocess_data.pk()] to each of its elements
#' @author Caroline Ring
#' @importFrom magrittr `%>%`
#' @export
#' @family methods for pk_faceted objects
preprocess_data.pk_faceted <- function(obj){

obj <- obj %>%
  dplyr::transmute(pk_object = purrr::map(pk_object,
                                          preprocess_data)
                   )

class(obj) <- c("pk_faceted", "pk", class(obj))
    return(obj)

}

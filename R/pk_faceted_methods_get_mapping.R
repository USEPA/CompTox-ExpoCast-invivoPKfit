#'get_mapping
#'
#'get_mapping for pk_faceted objects
#'
#'@param obj An object of class `pk_faceted`
#'@return A [tibble::tibble()] grouped by the faceting variables with a new variable ` get_mapping ` containing the results of applying [get_mapping.pk()] to each item in the list column of [pk()] objects, `obj$pk_object`
#'@export
#'@author Caroline Ring
#'@family methods for pk_faceted objects
get_mapping.pk_faceted <- function(obj, ...){
obj %>% dplyr::mutate(get_mapping =purrr::map(pk_object,get_mapping))
}

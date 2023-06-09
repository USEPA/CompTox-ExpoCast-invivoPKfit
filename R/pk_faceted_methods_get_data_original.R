#'get_data_original
#'
#'get_data_original for pk_faceted objects
#'
#'@param obj An object of class `pk_faceted`
#'@return A [tibble::tibble()] grouped by the faceting variables with a new variable ` get_data_original ` containing the results of applying [get_data_original.pk()] to each item in the list column of [pk()] objects, `obj$pk_object`
#'@export
#'@author Caroline Ring
#'@family methods for pk_faceted objects
get_data_original.pk_faceted <- function(obj, ...){
obj %>% dplyr::mutate(get_data_original =purrr::map(pk_object,get_data_original))
}

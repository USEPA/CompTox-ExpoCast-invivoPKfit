#'get_status
#'
#'get_status for pk_faceted objects
#'
#'@param obj An object of class `pk_faceted`
#'@return A [tibble::tibble()] grouped by the faceting variables with a new variable ` get_status ` containing the results of applying [get_status.pk()] to each item in the list column of [pk()] objects, `obj$pk_object`
#'@export
#'@author Caroline Ring
#'@family methods for pk_faceted objects
get_status.pk <- function(obj, ...){
obj %>% dplyr::mutate(get_status =purrr::map(pk_object,get_status))
}

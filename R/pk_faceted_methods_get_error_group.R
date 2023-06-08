#'get_error_group
#'
#'get_error_group for pk_faceted objects
#'
#'@param obj An object of class `pk_faceted`
#'@return A [tibble::tibble()] grouped by the faceting variables with a new variable ` get_error_group ` containing the results of applying [get_error_group.pk()] to each item in the list column of [pk()] objects, `obj$pk_object`
#'@export
#'@author Caroline Ring
#'@family methods for pk_faceted objects
get_error_group.pk <- function(obj, ...){
obj %>% dplyr::mutate(get_error_group =purrr::map(pk_object,get_error_group))
}

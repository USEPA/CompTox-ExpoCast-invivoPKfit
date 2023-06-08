#'get_data_sigma_group
#'
#'get_data_sigma_group for pk_faceted objects
#'
#'@param obj An object of class `pk_faceted`
#'@return A [tibble::tibble()] grouped by the faceting variables with a new variable ` get_data_sigma_group ` containing the results of applying [get_data_sigma_group.pk()] to each item in the list column of [pk()] objects, `obj$pk_object`
#'@export
#'@author Caroline Ring
#'@family methods for pk_faceted objects
get_data_sigma_group.pk <- function(obj, ...){
obj %>% dplyr::mutate(get_data_sigma_group =purrr::map(pk_object,get_data_sigma_group))
}

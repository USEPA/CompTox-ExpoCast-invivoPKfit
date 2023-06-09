#'get_settings_optimx
#'
#'get_settings_optimx for pk_faceted objects
#'
#'@param obj An object of class `pk_faceted`
#'@return A [tibble::tibble()] grouped by the faceting variables with a new variable ` get_settings_optimx ` containing the results of applying [get_settings_optimx.pk()] to each item in the list column of [pk()] objects, `obj$pk_object`
#'@export
#'@author Caroline Ring
#'@family methods for pk_faceted objects
get_settings_optimx.pk_faceted <- function(obj, ...){
obj %>% dplyr::mutate(get_settings_optimx =purrr::map(pk_object,get_settings_optimx))
}

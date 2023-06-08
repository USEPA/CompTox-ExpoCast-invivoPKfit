#'get_settings_preprocess
#'
#'get_settings_preprocess for pk_faceted objects
#'
#'@param obj An object of class `pk_faceted`
#'@return A [tibble::tibble()] grouped by the faceting variables with a new variable ` get_settings_preprocess ` containing the results of applying [get_settings_preprocess.pk()] to each item in the list column of [pk()] objects, `obj$pk_object`
#'@export
#'@author Caroline Ring
#'@family methods for pk_faceted objects
get_settings_preprocess.pk <- function(obj, ...){
obj %>% dplyr::mutate(get_settings_preprocess =purrr::map(pk_object,get_settings_preprocess))
}

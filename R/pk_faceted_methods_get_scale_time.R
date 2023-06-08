#'get_scale_time
#'
#'get_scale_time for pk_faceted objects
#'
#'@param obj An object of class `pk_faceted`
#'@return A [tibble::tibble()] grouped by the faceting variables with a new variable ` get_scale_time ` containing the results of applying [get_scale_time.pk()] to each item in the list column of [pk()] objects, `obj$pk_object`
#'@export
#'@author Caroline Ring
#'@family methods for pk_faceted objects
get_scale_time.pk <- function(obj, ...){
obj %>% dplyr::mutate(get_scale_time =purrr::map(pk_object,get_scale_time))
}

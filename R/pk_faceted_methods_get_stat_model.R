#'get_stat_model
#'
#'get_stat_model for pk_faceted objects
#'
#'@param obj An object of class `pk_faceted`
#'@return A [tibble::tibble()] grouped by the faceting variables with a new variable ` get_stat_model ` containing the results of applying [get_stat_model.pk()] to each item in the list column of [pk()] objects, `obj$pk_object`
#'@export
#'@author Caroline Ring
#'@family methods for pk_faceted objects
get_stat_model.pk <- function(obj, ...){
obj %>% dplyr::mutate(get_stat_model =purrr::map(pk_object,get_stat_model))
}

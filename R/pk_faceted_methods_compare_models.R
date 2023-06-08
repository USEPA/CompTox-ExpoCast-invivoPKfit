#'compare_models
#'
#'compare_models for pk_faceted objects
#'
#'@param obj An object of class `pk_faceted`
#'@return A [tibble::tibble()] grouped by the faceting variables with a new variable ` compare_models ` containing the results of applying [compare_models.pk()] to each item in the list column of [pk()] objects, `obj$pk_object`
#'@export
#'@author Caroline Ring
#'@family methods for pk_faceted objects
compare_models.pk <- function(obj, ...){
obj %>% dplyr::mutate(compare_models =purrr::map(pk_object,compare_models))
}

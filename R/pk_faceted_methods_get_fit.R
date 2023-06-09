#'get_fit
#'
#'get_fit for pk_faceted objects
#'
#'@param obj An object of class `pk_faceted`
#'@return A [tibble::tibble()] grouped by the faceting variables with a new variable ` get_fit ` containing the results of applying [get_fit.pk()] to each item in the list column of [pk()] objects, `obj$pk_object`
#'@export
#'@author Caroline Ring
#'@family methods for pk_faceted objects
get_fit.pk_faceted <- function(obj, ...){
obj %>% dplyr::mutate(get_fit =purrr::map(pk_object,get_fit))
}

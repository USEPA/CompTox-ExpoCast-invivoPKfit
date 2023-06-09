#'get_nca
#'
#'get_nca for pk_faceted objects
#'
#'@param obj An object of class `pk_faceted`
#'@return A [tibble::tibble()] grouped by the faceting variables with a new variable ` get_nca ` containing the results of applying [get_nca.pk()] to each item in the list column of [pk()] objects, `obj$pk_object`
#'@export
#'@author Caroline Ring
#'@family methods for pk_faceted objects
get_nca.pk_faceted <- function(obj, ...){
obj %>% dplyr::mutate(get_nca =purrr::map(pk_object,get_nca))
}

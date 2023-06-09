#'get_scale_conc
#'
#'get_scale_conc for pk_faceted objects
#'
#'@param obj An object of class `pk_faceted`
#'@return A [tibble::tibble()] grouped by the faceting variables with a new variable ` get_scale_conc ` containing the results of applying [get_scale_conc.pk()] to each item in the list column of [pk()] objects, `obj$pk_object`
#'@export
#'@author Caroline Ring
#'@family methods for pk_faceted objects
get_scale_conc.pk_faceted <- function(obj, ...){
obj %>% dplyr::mutate(get_scale_conc =purrr::map(pk_object,get_scale_conc))
}

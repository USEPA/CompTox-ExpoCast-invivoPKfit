#'subtract_pkproto
#'
#'subtract_pkproto for pk_faceted objects
#'
#'@param obj An object of class `pk_faceted`
#'@return A [tibble::tibble()] grouped by the faceting variables with a new variable ` subtract_pkproto ` containing the results of applying [subtract_pkproto.pk()] to each item in the list column of [pk()] objects, `obj$pk_object`
#'@export
#'@author Caroline Ring
#'@family methods for pk_faceted objects
subtract_pkproto.pk <- function(obj, ...){
obj %>% dplyr::mutate(subtract_pkproto =purrr::map(pk_object,subtract_pkproto))
}

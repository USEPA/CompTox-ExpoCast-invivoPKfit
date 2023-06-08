#'add_pkproto
#'
#'add_pkproto for pk_faceted objects
#'
#'@param obj An object of class `pk_faceted`
#'@return A [tibble::tibble()] grouped by the faceting variables with a new variable ` add_pkproto ` containing the results of applying [add_pkproto.pk()] to each item in the list column of [pk()] objects, `obj$pk_object`
#'@export
#'@author Caroline Ring
#'@family methods for pk_faceted objects
add_pkproto.pk <- function(obj, ...){
obj %>% dplyr::mutate(add_pkproto =purrr::map(pk_object,add_pkproto))
}

#'nca
#'
#'nca for pk_faceted objects
#'
#'@param obj An object of class `pk_faceted`
#'@return A [tibble::tibble()] grouped by the faceting variables with a new variable ` nca ` containing the results of applying [nca.pk()] to each item in the list column of [pk()] objects, `obj$pk_object`
#'@export
#'@author Caroline Ring
#'@family methods for pk_faceted objects
nca.pk_faceted <- function(obj, ...){
obj %>% dplyr::mutate(nca =purrr::map(pk_object,nca))
}

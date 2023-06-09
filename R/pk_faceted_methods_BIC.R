#'BIC
#'
#'BIC for pk_faceted objects
#'
#'@param obj An object of class `pk_faceted`
#'@return A [tibble::tibble()] grouped by the faceting variables with a new variable ` BIC ` containing the results of applying [BIC.pk()] to each item in the list column of [pk()] objects, `obj$pk_object`
#'@export
#'@author Caroline Ring
#'@family methods for pk_faceted objects
BIC.pk_faceted <- function(obj, ...){
obj %>% dplyr::mutate(BIC =purrr::map(pk_object,BIC))
}

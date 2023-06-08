#'coef
#'
#'coef for pk_faceted objects
#'
#'@param obj An object of class `pk_faceted`
#'@return A [tibble::tibble()] grouped by the faceting variables with a new variable ` coef ` containing the results of applying [coef.pk()] to each item in the list column of [pk()] objects, `obj$pk_object`
#'@export
#'@author Caroline Ring
#'@family methods for pk_faceted objects
coef.pk <- function(obj, ...){
obj %>% dplyr::mutate(coef =purrr::map(pk_object,coef))
}

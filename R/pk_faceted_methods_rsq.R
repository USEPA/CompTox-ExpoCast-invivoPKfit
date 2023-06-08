#'rsq
#'
#'rsq for pk_faceted objects
#'
#'@param obj An object of class `pk_faceted`
#'@return A [tibble::tibble()] grouped by the faceting variables with a new variable ` rsq ` containing the results of applying [rsq.pk()] to each item in the list column of [pk()] objects, `obj$pk_object`
#'@export
#'@author Caroline Ring
#'@family methods for pk_faceted objects
rsq.pk <- function(obj, ...){
obj %>% dplyr::mutate(rsq =purrr::map(pk_object,rsq))
}

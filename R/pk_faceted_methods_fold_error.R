#'fold_error
#'
#'fold_error for pk_faceted objects
#'
#'@param obj An object of class `pk_faceted`
#'@return A [tibble::tibble()] grouped by the faceting variables with a new variable ` fold_error ` containing the results of applying [fold_error.pk()] to each item in the list column of [pk()] objects, `obj$pk_object`
#'@export
#'@author Caroline Ring
#'@family methods for pk_faceted objects
fold_error.pk <- function(obj, ...){
obj %>% dplyr::mutate(fold_error =purrr::map(pk_object,fold_error))
}

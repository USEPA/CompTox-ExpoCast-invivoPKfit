#'predict
#'
#'predict for pk_faceted objects
#'
#'@param obj An object of class `pk_faceted`
#'@return A [tibble::tibble()] grouped by the faceting variables with a new variable ` predict ` containing the results of applying [predict.pk()] to each item in the list column of [pk()] objects, `obj$pk_object`
#'@export
#'@author Caroline Ring
#'@family methods for pk_faceted objects
predict.pk <- function(obj, ...){
obj %>% dplyr::mutate(predict =purrr::map(pk_object,predict))
}

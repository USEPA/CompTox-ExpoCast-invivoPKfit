#'rmse
#'
#'rmse for pk_faceted objects
#'
#'@param obj An object of class `pk_faceted`
#'@return A [tibble::tibble()] grouped by the faceting variables with a new variable ` rmse ` containing the results of applying [rmse.pk()] to each item in the list column of [pk()] objects, `obj$pk_object`
#'@export
#'@author Caroline Ring
#'@family methods for pk_faceted objects
rmse.pk_faceted <- function(obj, ...){
obj %>% dplyr::mutate(rmse =purrr::map(pk_object,rmse))
}

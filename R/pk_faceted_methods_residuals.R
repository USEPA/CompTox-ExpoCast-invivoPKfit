#'residuals
#'
#'residuals for pk_faceted objects
#'
#'@param obj An object of class `pk_faceted`
#'@return A [tibble::tibble()] grouped by the faceting variables with a new variable ` residuals ` containing the results of applying [residuals.pk()] to each item in the list column of [pk()] objects, `obj$pk_object`
#'@export
#'@author Caroline Ring
#'@family methods for pk_faceted objects
residuals.pk <- function(obj, ...){
obj %>% dplyr::mutate(residuals =purrr::map(pk_object,residuals))
}

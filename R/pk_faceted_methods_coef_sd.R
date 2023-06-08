#'coef_sd
#'
#'coef_sd for pk_faceted objects
#'
#'@param obj An object of class `pk_faceted`
#'@return A [tibble::tibble()] grouped by the faceting variables with a new variable ` coef_sd ` containing the results of applying [coef_sd.pk()] to each item in the list column of [pk()] objects, `obj$pk_object`
#'@export
#'@author Caroline Ring
#'@family methods for pk_faceted objects
coef_sd.pk <- function(obj, ...){
obj %>% dplyr::mutate(coef_sd =purrr::map(pk_object,coef_sd))
}

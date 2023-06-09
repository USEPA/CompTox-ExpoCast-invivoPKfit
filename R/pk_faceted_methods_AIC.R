#'AIC
#'
#'AIC for pk_faceted objects
#'
#'@param obj An object of class `pk_faceted`
#'@return A [tibble::tibble()] grouped by the faceting variables with a new variable ` AIC ` containing the results of applying [AIC.pk()] to each item in the list column of [pk()] objects, `obj$pk_object`
#'@export
#'@author Caroline Ring
#'@family methods for pk_faceted objects
AIC.pk_faceted <- function(obj, ...){
obj %>% dplyr::mutate(AIC =purrr::map(pk_object,AIC))
}

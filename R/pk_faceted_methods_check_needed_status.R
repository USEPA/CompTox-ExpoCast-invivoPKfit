#'check_needed_status
#'
#'check_needed_status for pk_faceted objects
#'
#'@param obj An object of class `pk_faceted`
#'@return A [tibble::tibble()] grouped by the faceting variables with a new variable ` check_needed_status ` containing the results of applying [check_needed_status.pk()] to each item in the list column of [pk()] objects, `obj$pk_object`
#'@export
#'@author Caroline Ring
#'@family methods for pk_faceted objects
check_needed_status.pk_faceted <- function(obj, ...){
obj %>% dplyr::mutate(check_needed_status =purrr::map(pk_object,check_needed_status))
}

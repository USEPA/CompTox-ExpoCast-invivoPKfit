#'data_summary
#'
#'data_summary for pk_faceted objects
#'
#'@param obj An object of class `pk_faceted`
#'@return A [tibble::tibble()] grouped by the faceting variables with a new variable ` data_summary ` containing the results of applying [data_summary.pk()] to each item in the list column of [pk()] objects, `obj$pk_object`
#'@export
#'@author Caroline Ring
#'@family methods for pk_faceted objects
data_summary.pk <- function(obj, ...){
obj %>% dplyr::mutate(data_summary =purrr::map(pk_object,data_summary))
}

#'plot
#'
#'plot for pk_faceted objects
#'
#'@param obj An object of class `pk_faceted`
#'@return A [tibble::tibble()] grouped by the faceting variables with a new variable ` plot ` containing the results of applying [plot.pk()] to each item in the list column of [pk()] objects, `obj$pk_object`
#'@export
#'@author Caroline Ring
#'@family methods for pk_faceted objects
plot.pk <- function(obj, ...){
obj %>% dplyr::mutate(plot =purrr::map(pk_object,plot))
}

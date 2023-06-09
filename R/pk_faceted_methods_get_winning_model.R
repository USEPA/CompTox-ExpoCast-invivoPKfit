#'get_winning_model
#'
#'get_winning_model for pk_faceted objects
#'
#'@param obj An object of class `pk_faceted`
#'@return A [tibble::tibble()] grouped by the faceting variables with a new variable ` get_winning_model ` containing the results of applying [get_winning_model.pk()] to each item in the list column of [pk()] objects, `obj$pk_object`
#'@export
#'@author Caroline Ring
#'@family methods for pk_faceted objects
get_winning_model.pk_faceted <- function(obj, ...){
obj %>% dplyr::mutate(get_winning_model =purrr::map(pk_object,get_winning_model))
}

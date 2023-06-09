#'eval_tkstats
#'
#'eval_tkstats for pk_faceted objects
#'
#'@param obj An object of class `pk_faceted`
#'@return A [tibble::tibble()] grouped by the faceting variables with a new variable ` eval_tkstats ` containing the results of applying [eval_tkstats.pk()] to each item in the list column of [pk()] objects, `obj$pk_object`
#'@export
#'@author Caroline Ring
#'@family methods for pk_faceted objects
eval_tkstats.pk_faceted <- function(obj, ...){
obj %>% dplyr::mutate(eval_tkstats =purrr::map(pk_object,eval_tkstats))
}

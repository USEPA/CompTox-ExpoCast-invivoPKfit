#' Get data_info
#'
#' Extract summary data information results from a [pk()] object
#'
#' @param obj A [pk()] object that has had `data_info()` run on it
#' @return A `data.frame`: the `data` element of `obj`
#' @export
#' @author Caroline Ring
get_data_info.pk <- function(obj){
  #check if data has been data_info(0)
  check <- check_required_status(obj = obj,
                                 required_status = status_data_info)
  if(!(check %in% TRUE)){
    stop(attr(check, "msg"))
  }

  return(obj$data_info)
}

#' Get NCA
#'
#' Extract non-compartmental analysis results from a [pk()] object
#'
#' @param obj A [pk()] object that has had `data_info()` run on it
#' @return A `data.frame`: the `data` element of `obj`
#' @export
#' @author Caroline Ring, Gilberto Padilla Mercado
get_nca.pk <- function(obj){
  #check if data has been data_info(0)
  check <- check_required_status(obj = obj,
                                 required_status = status_data_info)
  if(!(check %in% TRUE)){
    stop(attr(check, "msg"))
  }

  return(obj$data_info$nca)
}

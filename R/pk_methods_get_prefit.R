#' Get prefit
#'
#' Extract pre-fitting results from a [pk()] object
#'
#' @param obj A [pk()] object that has had `data_info()` run on it
#' @return A list of `data.frame`s: one for each model in `obj$stat_model`
#' @export
#' @author Caroline Ring
get_prefit.pk <- function(obj,
                          model = NULL){

  #check if data has been prefit
  check <- check_required_status(obj = obj,
                                 required_status = status_prefit)
  if(!(check %in% TRUE)){
    stop(attr(check, "msg"))
  }

  if(is.null(model)){
    model <- names(obj$stat_model)
  }

  parDF_list <- sapply(model,
                function(this_model) {
                  tmp <- rbind(
                    obj$prefit[[this_model]]$par_DF,
                    obj$prefit$stat_error_model$sigma_DF)
                  rownames(tmp) <- NULL
                  tmp
                  },
                simplify = FALSE,
                USE.NAMES = TRUE)

  return(parDF_list)
}

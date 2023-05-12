#' Plot a [pk()] object.
#'
#' Plot data and model fits from a [pk()] object.
#'
#' If the [pk()] object has not been fitted, then only the data will be plotted
#' (because no curve fits exist).
#'
#' @param obj A [pk()] object
#' @param newdata Optional: A `data.frame` containing new data to plot. Must
#'   contain at least variables `Chemical`, `Species`, `Route`, `Media`, `Dose`,
#'   `Time`, `Time.Units`, `Conc`, `Detect`, `Conc_SD`.  Default `NULL`, to use
#'   the data in `obj$data`.
#' @param model Character: One or more of the models fitted. Curve fits will be
#'   plotted for these models. Default `NULL` to plot fits for all models in
#'   `obj$stat_model`.
#' @param method Character: One or more of the [optimx::optimx()] methods used.
#'   Default `NULL` to plot fits for all methods in
#'   `obj$optimx_settings$method`.
plot.pk <- function(obj,
                    newdata = NULL,
                    model = NULL,
                    method = NULL){

}

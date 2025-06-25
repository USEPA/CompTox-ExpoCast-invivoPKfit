#' Get stat_model
#'
#' @param obj A [pk()] object
#' @param ... Additional arguments. Currently not in use.
#' @return A `list` -- the `stat_model` element of `obj`
#' @export
#' @author Caroline Ring
get_stat_model.pk <- function(obj, ...) {

  tidyr::tibble(modelfun = obj$stat_model) |>
    dplyr::mutate(model = purrr::map(modelfun, \(x) {x$name})) |>
    tidyr::unnest(model)

}

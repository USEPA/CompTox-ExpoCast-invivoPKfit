#'Plot a [pk()] object.
#'
#'Plot data and model fits from a [pk()] object.
#'
#'If the [pk()] object has not been fitted, then only the data will be plotted
#'(because no curve fits exist).
#'
#'@param obj A [pk()] object
#'@param newdata Optional: A `data.frame` containing new data to plot. Must
#'  contain at least variables `Chemical`, `Species`, `Route`, `Media`, `Dose`,
#'  `Time`, `Time.Units`, `Conc`, `Detect`, `Conc_SD`.  Default `NULL`, to use
#'  the data in `obj$data`.
#'@param model Character: One or more of the models fitted. Curve fits will be
#'  plotted for these models. Default `NULL` to plot fits for all models in
#'  `obj$stat_model`.
#'@param method Character: One or more of the [optimx::optimx()] methods used.
#'  Default `NULL` to plot fits for all methods in `obj$settings_optimx$method`.
#'@param use_scale_conc Possible values: `TRUE`, `FALSE`, or a named list with
#'  elements `dose_norm` and `log10_trans` which themselves should be either
#'  `TRUE` or `FALSE`.  If `use_scale_conc = FALSE` (the default for this
#'  function), then the data and fits will be plotted without any
#'  dose-normalization or log-transformation. If `use_scale_conc = TRUE` , then
#'  the concentration scaling/transformations in `obj` will be applied to the
#'  y-axis (concentration axis). If `use_scale_conc = list(dose_norm = ...,
#'  log10_trans = ...)`, then the specified dose normalization and/or
#'  log10-transformation will be applied to the y-axis (concentration axis) of
#'  the plots.
#'@param plot_data_aes Optional: Aesthetic mapping for the plot layer that
#'  visualizes the data. Default `NULL`, in which case a default mapping will be
#'  used based on the value of `use_scale_conc`.
#'@param plot_fit_aes Optional: Aesthetic mapping for the plot layer that
#'  visualizes the fitted curves. Default `NULL`, in which case a default mapping will be
#'  used based on the value of `use_scale_conc`.
#'@param facet_fun  Default `"facet_grid"`. Optional: The name of the `ggplot2` faceting function to use:
#'  [ggplot2::facet_grid()], [ggplot2::facet_wrap()], or `'none'` to do no
#'  faceting. Default `NULL`, in which case a default faceting will be
#'  applied based on the value of `use_scale_conc`.
#'@param facet_fun_args A named list of arguments to the faceting function in
#'  `facet_fun` (if any). Default:
#' ```
#'list(rows = ggplot2::vars(Route),
#'cols= ggplot2::vars(Media),
#'scales = "free_y",
#'labeller = "label_both")
#'```
#'@param n_interp For plotting: the number of time points to interpolate between
#'  each observed time point. Default 10.
#'@return A [ggplot2::ggplot()]-class plot object.
#'@import ggplot2
#'@export
plot.pk <- function(obj,
                    newdata = NULL,
                    model = NULL,
                    method = NULL,
                    use_scale_conc = FALSE,
                    # Plotting arguments
                    plot_data_aes = NULL,
                    plot_point_aes = NULL,
                    facet_fun = NULL,
                    facet_fun_args = NULL,
                    drop_nonDetect = FALSE,
                    # Predict/interpolation arguments
                    plot_fit_aes = NULL,
                    n_interp = 10,
                    ...){

  #ensure that the model has at least been preprocessed
  check <- check_required_status(obj = obj,
                                 required_status = status_preprocess)
  if (!(check %in% TRUE)){
    stop(attr(check, "msg"))
  }

  if (is.null(model)) model <- names(obj$stat_model)
  if (is.null(method)) method <- obj$settings_optimx$method
  if (is.null(newdata)) newdata <- obj$data

  newdata_ok <- check_newdata(newdata = newdata,
                              olddata = obj$data,
                              req_vars = c("Time",
                                           "Time.Units",
                                           "Dose",
                                           "Route",
                                           "Media"),
                              exclude = FALSE)

  #check method and model
  method_ok <- check_method(obj = obj, method = method)
  model_ok <- check_model(obj = obj, model = model)

  #apply transformations if so specified
  conc_scale <- conc_scale_use(obj = obj,
                               use_scale_conc = use_scale_conc)

  if (drop_nonDetect %in% TRUE) {
    newdata <- subset(newdata, Detect %in% TRUE)
  }

  req_vars <- c("Time",
                "Time.Units",
                "Dose",
                "Route",
                "Media",
                "Reference",
                "Conc",
                "Conc_SD",
                "Value",
                "Value.Units",
                "Detect", "exclude")

  trans_vars <- c("Time_trans",
                  "Time_trans.Units",
                  "Conc_trans",
                  "Conc_trans.Units")

  # I think ideally there should be required variables and
  # a way to ensure all the aes() variables get added
  newdata <- newdata %>%
    dplyr::select(!!!obj$data_group,
                  dplyr::all_of(req_vars),
                  dplyr::any_of(trans_vars)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(Conc_set = ifelse(conc_scale$dose_norm,
                                    Conc / Dose,
                                    Conc),
                  Conc_set_SD = ifelse(conc_scale$dose_norm,
                                       Conc_SD / Dose,
                                       Conc_SD),
                  Detect = ifelse(Detect,
                                  "Detect",
                                  "Non-Detect"),
                  Dose = factor(Dose)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(!!!obj$data_group) %>%
    tidyr::nest(.key = "observations")

  # Default arguments and functions
  plot_data_aes_default <- ggplot2::aes(x = Time_trans,
                                        y = Conc_set,
                                        ymin = Conc_set - Conc_set_SD,
                                        ymax = Conc_set + Conc_set_SD,
                                        color = Dose,
                                        shape = Reference)

  plot_point_aes_default <- ggplot2::aes(fill = Dose,
                                         alpha = Detect)

  plot_fit_aes_default <- ggplot2::aes(x = Time_trans,
                                       y = Conc_est,
                                       linetype = interaction(model, method))

  facet_fun_default <- ifelse(conc_scale$dose_norm,
                              "facet_grid", "facet_wrap")

  facet_fun_args_default <- ifelse(conc_scale$dose_norm,
                                   list(rows = ggplot2::vars(Route),
                                        cols = ggplot2::vars(Media),
                                        scales = "free_y",
                                        labeller = "label_both"),
                                   list(ggplot2::vars(Route, Media, Dose),
                                        scales = "free_y",
                                        labeller = "label_both"))

  if (is.null(plot_data_aes)) plot_data_aes <- plot_data_aes_default
  if (is.null(plot_point_aes)) plot_point_aes <- plot_point_aes_default
  if (is.null(plot_fit_aes)) plot_fit_aes <- plot_fit_aes_default
  if (is.null(facet_fun)) facet_fun <- facet_fun_default
  if (is.null(facet_fun_args)) facet_fun_args <- facet_fun_args_default

  if (facet_fun %in% "none") facet_fun <- NULL

  newdata <- newdata %>%
    mutate(observation_plot =
             map(observations,
                 \(x) {
                   # Need to write this into a mutate + map pattern function
                   #initialize plot
                   p <- ggplot(data = x,
                               mapping = plot_data_aes) +
                     geom_point(mapping = plot_point_aes,
                                stroke = 1) + # plots points
                     geom_errorbar(mapping = plot_data_aes)

                   # if alpha = Detect then we have to do some trickery to implement that
                   # that is the objective of plot_point_aes
                   if (rlang::as_label(plot_point_aes$alpha) %in% "Detect") {
                     # plot the data
                     p <- p +
                       geom_point()
                     # plots here are un-filled points (only inherits color variable by default)

                     if ("shape" %in% c(names(plot_data_aes),
                                        names(plot_point_aes))) {
                       # Ensure shapes can take both fill and color if alpha is set
                       p <- p +
                         scale_shape_manual(values = 21:25)
                     }

                     #set an alpha scale
                     p <- p + scale_alpha_manual(values = c("Detect" = 1, "Non-Detect" = 0),
                                                 breaks = c("Detect", "Non-Detect"),
                                                 drop = FALSE,
                                                 name = NULL)

                     # Need to check whether there are two values in Detect
                     # some might be all Non-Detect or all Detect
                     if (length(unique(x$Detect)) > 1){
                       alpha_fill <- c("black", NA)
                     } else {
                       alpha_fill <- ifelse(unique(x$Detect) == "Detect",
                                            "black", NA)
                     }
                     p <- p +
                       guides(alpha = guide_legend(override.aes = list(shape = 21,
                                                                       color = "black",
                                                                       stroke = 1,
                                                                       fill = alpha_fill,
                                                                       alpha = 1)))

                   }

                   #Apply faceting, if any
                   if (!is.null(facet_fun)) {
                     p <- p +
                       do.call(facet_fun,
                               args = facet_fun_args)
                   }

                   p +
                     labs(title = paste(Chemical, Species),
                          subtitle = method,
                          x = "Time",
                          y = "Concentration") +
                     theme_bw() +
                     theme(panel.border = element_rect(color = "black", fill = NA,
                                                       linewidth = 1))

                 }))

  # For predictions, I will interpolate the time

  if (get_status(obj = obj) == 5) {
    interp_data <- newdata
    interp_data <- interp_data %>%
      mutate(interpolated = map(observations,
                                \(x) {
                                  x %>%
                                    dplyr::select(Time, Time.Units,
                                                  Dose, Route, Media) %>%
                                    dplyr::distinct() %>%
                                    tidyr::uncount(n_interp) %>%
                                    dplyr::group_by(Dose, Route, Media) %>%
                                    dplyr::mutate(Time_trans = (max(Time) / (n() - 1)) *
                                                    (row_number() - 1),
                                                  Dose = as.numeric(Dose),
                                                  Chemical = Chemical,
                                                  Species = Species)
                                }))
    interp_data <- interp_data %>%
      mutate(predicted = map(interpolated,
                             \(x) {
                               predict(obj = obj,
                                       newdata = x,
                                       use_scale_conc = conc_scale$dose_norm,
                                       exclude = FALSE)
                             }))


    return(interp_data)
  }

  return(newdata)
}

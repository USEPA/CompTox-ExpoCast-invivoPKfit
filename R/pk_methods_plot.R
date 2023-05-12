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
#' @param plot_data_aes Aesthetic mapping for the plot layer that visualizes the
#'   data.
#' @param plot_fit_aes Aesthetic mapping for the plot layer that visualizes the
#'   fits.
#' @param facet_fun The name of the `ggplot2` faceting function to use:
#'   [ggplot2::facet_grid()], [ggplot2::facet_wrap()], or `NULL` to do no
#'   faceting. Default `"facet_grid"`.
#' @param facet_fun_args A named list of arguments to the faceting function in `facet_fun`. Default:
#' ```
#'list(rows = ggplot2::vars(Route),
#'cols= ggplot2::vars(Media),
#'scales = "free_y",
#'labeller = "label_both")
#'```
#' @param n_interp For plotting: the number of time points to interpolate between each observed time point. Default 10.
#' @return A [ggplot2::ggplot()]-class plot object.
#' @import ggplot2
#' @export
plot.pk <- function(obj,
                    newdata = NULL,
                    model = NULL,
                    method = NULL,
                    plot_data_aes = ggplot2::aes(x = Time_trans,
                                                y = Conc/Dose,
                                                color = as.factor(Dose),
                                                ymin = Conc/Dose - Conc_SD/Dose,
                                                ymax = Conc/Dose + Conc_SD/Dose,
                                                shape = Reference,
                                                fill = Detect),
                    plot_fit_aes = ggplot2::aes(x = Time_trans,
                                                y = Conc/Dose,
                                                linetype = interaction(model, method)),
                    facet_fun = "facet_grid",
                    facet_fun_args = list(rows = ggplot2::vars(Route),
                                           cols= ggplot2::vars(Media),
                                           scales = "free_y",
                                           labeller = "label_both"),
                    n_interp = 10
                    ){

  #ensure that the model has been fitted
  check <- check_required_status(obj = obj,
                                 required_status = status_fit)
  if(!(check %in% TRUE)){
    stop(attr(check, "msg"))
  }

  if(is.null(model)) model <- names(obj$stat_model)
  if(is.null(method)) method <- obj$optimx_settings$method
  if(is.null(newdata)) newdata <- obj$data

  #check that all methods are valid
  if(!(all(method %in% obj$optimx_settings$method))){
    stop(paste("All values in `method` must be found in `obj$optimx_settings$method.",
               paste0("`method` = ", paste(method, sep = ", ")),
               paste0("`obj$optimx_settings$method` = ", paste(obj$optimx_settings$method)),
               sep = "\n"))
  }

  p <- ggplot(data = newdata) +
    geom_blank() +
    theme_bw()

  #first plot data
  if (rlang::as_label(plot_data_aes$fill) %in% "Detect"){
    newdata$Detect <- factor(newdata$Detect,
                             levels = c(TRUE, FALSE),
                             labels = c("Detect", "Non-Detect"))
    #if fill = Detect then we have to do some trickery to implement that
    #make sure there is a shape aesthetic -- create one if there is no
    if(!("shape" %in% names(plot_data_aes))){
      no_shape <- TRUE
      plot_data_aes <- c(plot_data_aes,
                         aes(shape = "junk"))
    }else{
      no_shape <- FALSE
    }
    #replace "fill" aesthetic with "colour" aesthetic, if any
    #have to use spelling "colour" with a U, because that's how it's spelled in the aes object list names
    plot_data_aes$fill <- plot_data_aes$colour
    #create a mapping with alpha mapped to detect
    #do this by concatenating two aes() mappings
    alpha_map <- c(plot_data_aes,
                   ggplot2::aes(alpha = Detect))
    #the above produces a list; ggplot2 doesn't like it beacuse it doesn't have class "uneval", even though it is structurally the same thing
    #so, give it that class
    class(alpha_map) <- c(class(alpha_map), "uneval") #the class expected by ggplot2
    #make a tmp aesthetic without any fill
    aes_tmp <- copy(plot_data_aes)
    aes_tmp$fill <- NULL

    if("shape" %in% names(plot_data_aes)){
      p <- p +
        scale_shape_manual(values = 21:25)
    }

    #set an alpha scale
    p <- p + scale_alpha_manual(values = c("Detect" = 1, "Non-Detect" = 0),
                                drop = FALSE,
                                name = NULL)

    p <- p +
      geom_errorbar(mapping = plot_data_aes)


    p <- p  +
      geom_point(data = subset(newdata,
                               Detect %in% "Detect"),
                 mapping = aes_tmp, #first plot detects with white fill
                 fill = "white",
                 size = 4,
                 stroke = 1)
    p <- p +
      geom_point(data = subset(newdata,
                               Detect %in% "Detect"),
                 mapping = alpha_map, #then plot detects with fill, and alpha mapped to Detect
                 size = 4,
                 stroke = 1)

    p <- p +
      geom_point(data = subset(newdata,
                               Detect %in% "Non-Detect"), #then plot non-detects with white fill
                 mapping = aes_tmp,
                 fill = "white",
                 size = 4,
                 stroke = 1)

    p <- p +
      geom_point(data = subset(newdata,
                               Detect %in% "Non-Detect"),
                 mapping = alpha_map, #then plot non-detects with fill, and alpha mapped to Detect
                 size = 4,
                 stroke = 1)


p <- p +
      guides(alpha = guide_legend(override.aes = list(shape = 21,
                                                      color = "black",
                                                      stroke = 1,
                                                      fill = c("black", NA),
                                                      alpha = 1)))
    if(no_shape){
      #suppress legend for shape if there was not originally any shape aesthetic
      p <- p + guides(shape = "none")
    }
  }else{
    p <- p +
      geom_errorbar(mapping = plot_data_aes) +
      geom_point(mapping = plot_data_aes,
                 stroke = 1)
  }


p <- p +
    do.call(facet_fun,
            args = facet_fun_args)



  #then create a data.frame of predicted concentrations at interpolated time points

  data_plot <- newdata %>%
    dplyr::group_by(Chemical, Species, Route, Media, Dose, Reference) %>%
    dplyr::summarise(Time_trans = {
      #get unique timepoints in this group, including 0, sorted in increasing order
      time_tmp <- sort(unique(c(0,Time_trans)))
      #between each of those time points, interpolate 10 new points
      time_interp <- sapply(1:(length(time_tmp)-1),
                            function(i) {
                              seq(from = time_tmp[i],
                                  to = time_tmp[i+1],
                                  length.out = n_interp)
                            })
      #convert from matrix into vector, and keep only unique items
      unique(as.numeric(time_interp))
    }
    ) %>%
    dplyr::ungroup() %>%
    as.data.frame()

  #for each model, get predictions for the plot data
  preds <- predict(obj = obj,
                   newdata = data_plot,
                   model = model,
                   method = method,
                   type = "conc")
  #translate the predictions into data.frames with variables model, method, and Conc
  preds_all <- do.call(rbind,
          lapply(model,
                 function(this_model){
                   do.call(rbind,
                           lapply(method,
                                  function(this_method){
                                    preds_vect <- preds[[this_model]][, this_method]
                                    preds_DF <- data.frame(model = this_model,
                                               method = this_method,
                                               Conc = preds_vect)
                                    #column bind these predictions to data_plot,
                                    #and return the result
                                    cbind(data_plot,
                                          preds_DF)
                                  }
                           )
                   )
                 }
          )
  )

  #now, add the predictions layer
  p <- p +
    geom_line(data = preds_all,
              mapping = plot_fit_aes)

  return(p)
}

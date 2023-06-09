#' Add a `pkproto` object to a `pk` object
#'
#' Add a `pkproto` object to a `pk` object
#'
#' Note that `e1 + e2` is equivalent to
#'
#' ```
#' `+`(e1, e2)
#' ```
#'
#' @param e1 A `pk` pbject
#' @param e2 A `pkproto` object
#' @return The `pk` object, modified by adding the `pkproto` object
#'@export
#' @author Caroline Ring
"+.pk" <- function(e1, e2) {
  #throw error if only one argument
  if (missing(e2)) {
    cli::cli_abort(c(
      "Cannot use {.code +} with a single argument",
      "i" = "Did you accidentally put {.code +} on a new line?"
    ))
  }

  # Get the name of what was passed in as e2, and pass along so that it
  # can be displayed in error messages
  e2name <- deparse(substitute(e2))

  #make sue
  if      (is.pk(e1))  add_pk(e1, e2, e2name)
  else if (is.pkproto(e1)) {
    cli::cli_abort(c(
      "Cannot add {.cls pkproto} objects together",
      "i" = "Did you forget to add this object to a {.cls pk} object?"
    ))
  }
}

#' Add various `pkproto` objects to a `pk` object
#'
#'@param pk_obj The `pk` object
#'@param object The `pkproto` object to be added
#'@param objectname The name of the `pkproto` object to be added
#'
#'@return The `pk` object modified by the addition.
#' @export
add_pk <- function(pk_obj, object, objectname) {
  if (is.null(object)) return(pk_obj)

  p <- pk_add(object, pk_obj, objectname)
  p
}

#' Add a `pkproto` object to a `pk` object
#'
#' This is the S3 generic method.
#'
#' @param pkproto_obj The `pkproto` object to be added
#' @param pk_obj The `pk` object to which the `pkproto` object is to be added
#' @param objectname The object name
#'
#' @export
#' @seealso [pk_add.pk_scales()] for the method for adding `pk_scales` objects
#'   (from [scale_conc()] and [scale_time()]); [pk_add.pk_settings_preprocess()]
#'   for the method for adding `pk_settings_preprocess` objects (from
#'   [settings_preprocess()]); [pk_add.pk_settings_data_info()] for the method
#'   for adding `pk_settings_data_info` objects (from [settings_data_info()]);
#'   [pk_add.pk_settings_optimx()] for the method for adding
#'   `pk_settings_optimx` objects (from [settings_optimx()]);
#'   [pk_add.pk_stat_model()] for the method for adding `pk_stat_model` objects
#'   (from `stat_model()`)
pk_add <- function(pkproto_obj, pk_obj, objectname){
  UseMethod("pk_add")
}

#' Add pkproto object default method
#'
#' @export
pk_add.default <- function(object, pk_obj, objectname){
  stop(paste("No 'pk_add' method exists for object of class",
             object))
}

#'Subtract a `pkproto` object from a `pk` object
#'
#' This is the S3 generic method.
#'
#' @param pkproto_obj The `pkproto` object to be subtracted
#' @param pk_obj The `pk` object to which the `pkproto` object is to be subtracted
#' @param objectname The object name
#'
#' @export
#' @seealso
#'   [pk_subtract.pk_stat_model()] for the method for subtracting `pk_stat_model` objects
#'   (from `stat_model()`)
pk_subtract <- function(pkproto_obj, pk_obj, objectname){
  UseMethod("pk_subtract")
}

#' Subtract pkproto object default method
#'
#' @export
pk_subtract.default <- function(object, pk_obj, objectname){
  stop(paste("No 'pk_subtract' method exists for object of class",
             object))
}

#' Add a `pk_scales` object to a `pk` object.
#'
#' @param object The `pk_scales` object to be added.
#' @param pk_obj The `pk` object to which the `pk_scales` object will be added.
#' @param objectname The name of the `pk_scales` object.
#'
#' @return The `pk` object, modified by the `pk_scales` object.
#' @author Caroline Ring
#' @export
pk_add.pk_scales <- function(object, pk_obj, objectname){
  pk_obj$scales[[object$name]] <- object$value
  if(pk_obj$status > (status_preprocess - 1)){
    #with new scaling, everything will change starting from data pre-processing
    #with new data pre-processing settings, everything will change starting from
    #data pre-processing
    pk_obj$status <- status_preprocess - 1
    message(paste0(objectname,
                   ": Modifying data scaling resets status to level ",
                   status_preprocess - 1,
                   "; all later stages of the workflow will need to be re-done")
    )
    pk_obj$status <- status_preprocess - 1
  }

  return(pk_obj)
}

#' Add a `pk_settings_preprocess` object.
#'
#' @param object The `pk_settings_preprocess` object to be added.
#' @param pk_obj The `pk` object to which the `pk_settings_preprocess` object will be added.
#' @param objectname The name of the `pk_settings_preprocess` object.
#'
#' @return The `pk` object, modified by the `pk_settings_preprocess` object.
#' @author Caroline Ring
#' @export
pk_add.pk_settings_preprocess <- function(object, pk_obj, objectname){

  #New settings_preprocess will *replace* existing ones
  if(!is.null(pk_obj$settings_preprocess)){
    message(paste0(objectname,
                   ": settings_preprocess already present; new settings_preprocess will replace the existing one")
    )
  }

  pk_obj$settings_preprocess <- object
  if(pk_obj$status > (status_preprocess - 1)){
    #with new data pre-processing settings, everything will change starting from
    #data pre-processing
    pk_obj$status <- status_preprocess - 1
    message(paste0(objectname,
                   ": Modifying settings_preprocess resets status to level ",
                   status_preprocess - 1,
                   "; all later stages of the workflow will need to be re-done")
    )
  }
  return(pk_obj)
}

#' Add a `pk_settings_data_info` object.
#'
#' @param object The `pk_settings_data_info` object to be added.
#' @param pk_obj The `pk` object to which the `pk_settings_data_info` object will be added.
#' @param objectname The name of the `pk_settings_data_info` object.
#'
#' @return The `pk` object, modified by the `pk_settings_data_info` object.
#' @author Caroline Ring
#' @export
pk_add.pk_settings_data_info <- function(object, pk_obj, objectname){

  #New settings_data_info will *replace* existing ones
  if(!is.null(pk_obj$settings_data_info)){
    message(paste0(objectname,
                   ": settings_data_info already present; new settings_data_info will replace the existing one")
    )
  }

  pk_obj$settings_data_info <- object
  if(pk_obj$status > (status_preprocess - 1)){
    #with new data info settings, everything will change starting from
    #data info
    pk_obj$status <- status_data_info - 1
    message(paste0(objectname,
                   ": Modifying settings_data_info resets status to level ",
                   status_data_info - 1,
                   "; all later stages of the workflow will need to be re-done")
    )
  }
  return(pk_obj)
}

#' Add a `pk_settings_optimx` object.
#'
#' @param object The `pk_settings_optimx` object to be added.
#' @param pk_obj The `pk` object to which the `pk_settings_optimx` object will be added.
#' @param objectname The name of the `pk_settings_optimx` object.
#'
#' @return The `pk` object, modified by adding the settings.
#' @author Caroline Ring
#' @export
pk_add.pk_settings_optimx <- function(object, pk_obj, objectname){

  #New settings_optimx will *replace* existing ones
  if(!is.null(pk_obj$settings_optimx)){
    message(paste0(objectname,
                   ": settings_optimx already present; new settings_optimx will replace the existing one")
    )
  }

  pk_obj$settings_optimx <- object
  if(pk_obj$status > (status_prefit - 1)){
    #with new optimizer settings, data pre=processing and model pre-fitting
    #should not change, but model fitting will change
    message(paste0(objectname,
                   ": Modifying settings_optimx resets status to level ",
                   (status_prefit - 1),
                  "; all later stages of the workflow will need to be re-done")
    )
    pk_obj$status <- 3L
  }
  return(pk_obj)
}


#' Add a `pk_stat_model` object.
#'
#' @param object The `pk_stat_model` object to be added.
#' @param pk_obj The `pk` object to which the `pk_stat_model` object will be added.
#' @param objectname The name of the `pk_stat_model` object.
#'
#' @return The `pk` object, modified by adding the `stat_model`.
#' @author Caroline Ring
#' @export
pk_add.pk_stat_model <- function(object, pk_obj, objectname){
  #New stat_model will *replace* existing stat_model
  if(!is.null(pk_obj$stat_model)){
    old_models <- names(pk_obj$stat_model)
    new_models <- names(object)
    message(paste0(objectname,
                   ": stat_model already present; ",
                   "new stat_model will replace the existing one.",
                   "\n",
                   paste("Old models: ", paste(old_models, collapse = ", ")),
                   "\n",
                   paste("New models: ", paste(new_models, collapse = ", "))
                   ))
  }
  pk_obj$stat_model <- object

  #data pre-processing won't change with addition of new model, but model
  #pre-fit and fit will change
  if(pk_obj$status > (status_prefit - 1)){
    message(paste0(objectname,
                   ": Modifying stat_model resets status to level ",
                   status_prefit - 1,
                   "; all later stages of the workflow will need to be re-done")
    )
    pk_obj$status <- (status_prefit - 1)
  }

  return(pk_obj)
}

#' Add a `pk_stat_error_model` object.
#'
#' @param object The `pk_stat_error_model` object to be added.
#' @param pk_obj The `pk` object to which the `pk_stat_error_model` object will be added.
#' @param objectname The name of the `pk_stat_error_model` object.
#'
#' @return The `pk` object, modified by adding the `stat_error_model`.
#' @author Caroline Ring
#' @export
pk_add.pk_stat_error_model <- function(object, pk_obj, objectname){

  if(!is.null(pk_obj$stat_error_model)){
    message(paste0(objectname,
                   ": stat_error_model already present; new stat_error_model will replace the existing one")
    )
  }

  pk_obj$stat_error_model <- object

  #data pre-processing won't change with addition of a new error model, but
  #model pre-fit and fit will change
  if(pk_obj$status > (status_prefit - 1)){
    message(paste0(objectname,
                   ": Modifying stat_error_model resets status to level ",
                   status_prefit - 1,
                   "; all later stages of the workflow will need to be re-done"))
    pk_obj$status <-  status_prefit - 1
  }

  return(pk_obj)
}

#' Add an `uneval` object
#'
#' Add an object created by [ggplot2::aes()]
#'
#' This function adds a new variable mapping (created by [ggplot2::aes()]),
#' which has class `uneval`, to an existing [pk()] object.
#'
#' The new mapping will completely replace any existing mapping.
#'
#' @param object The `uneval` (mapping) object to be added.
#' @param pk_obj The [pk()] object to which the `uneval` object will be added.
#' @param objectname The name of the `uneval` object.
#'
#' @return The [pk()] object, modified by adding the new mapping.
#' @author Caroline Ring
#' @export
pk_add.uneval <- function(object, pk_obj, objectname){
  if(!is.null(pk_obj$mapping)){
    message(paste0(objectname,
                   ": mapping already present; new mapping will replace the existing one")
    )
  }

  pk_obj$mapping <- object

  #data pre-processing and everything downstream will change
  if(pk_obj$status > (status_preprocess - 1)){
    message(paste0(objectname,
                   ": Modifying mapping resets status to level ",
                   status_preprocess - 1,
                   "; all later stages of the workflow will need to be re-done"))
    pk_obj$status <-  status_preprocess - 1
  }

  return(pk_obj)
}

#' Add facet_data()
#'
#' @param object The `pk_facet_data` object to be added.
#' @param pk_obj The [pk()] object to which the `pk_facet_data` object will be added.
#' @param objectname The name of the `pk_facet_data` object.
pk_add.pk_facet_data <- function(object, pk_obj, objectname){

  #data pre-processing and everything downstream will change
  if(pk_obj$status > (status_preprocess - 1)){
    message(paste0(objectname,
                   ": Modifying faceting resets status to level ",
                   status_preprocess - 1,
                   "; all later stages of the workflow will need to be re-done"))
    pk_obj$status <-  status_preprocess - 1
  }

  #this will become the data_group
  mapping <- sapply(pk_obj$mapping, rlang::as_label)
  facets <- sapply(object$facets, rlang::as_label)
  newvars <- names(mapping)[match(facets,mapping)]
  #construct data_group as list of quosures
  #(same effect as a `dplyr::vars()` call)
  data_group <-  do.call(rlang::quos, rlang::syms(newvars))

  #check to see whether the facets produce more than one group of data
  #if they produce only one, then do not change class to pk_faceted

  data_original <- do.call(dplyr::group_by,
                      args = c(list(pk_obj$data_original),
                      object$facets))

  n_grps <- dplyr::n_groups(data_original)

  if(n_grps > 1){
    #translate the faceting variables (original data variables names) into harmonized variable names


    #produce a tibble with a list column of pk objects, each for one slice of the data
    pk_faceted_obj <-  do.call(dplyr::group_by,
                               args = c(list(pk_obj$data_original),
                                        object$facets)) %>%
      dplyr::summarise(pk_object = list(
        pk(data = dplyr::cur_data_all(), #the current slice of data
                                   mapping = pk_obj$mapping, #use the same mapping
           data_group = data_group,
           #use the same settings, scales, stats
           settings_preprocess_args = get_settings_preprocess(pk_obj),
           settings_data_info_args = get_settings_data_info(pk_obj),
           settings_optimx_args = get_settings_optimx(pk_obj),
           scale_conc_args = get_scale_conc(pk_obj),
           scale_time_args = get_scale_time(pk_obj),
           stat_model_args = get_stat_model(pk_obj),
           stat_error_model_args = get_stat_error_model(pk_obj))
                                   )
        )
#set an attribute to store faceting variables
    attr(pk_faceted_obj, "facets") <- object$facets
    class(pk_faceted_obj) <- c("pk_faceted",
                               "pk",
                               class(pk_faceted_obj))
    return(pk_faceted_obj)
  }else{
    #if there is only one facet,
    #just keep treating it as a `pk` object -- no need to go to all the trouble
    #just modify data_group to match faceting
    pk_obj$data_group <- data_group
    return(pk_obj)
  }
}

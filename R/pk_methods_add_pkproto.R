#' Add a pkproto object to a pk object
#'
#' @param e1 A pk pbject
#' @param e2 A pkproto object
#' @return The pk object, modified by adding the pkproto object
#'@export
#' @author Caroline Ring
"+.pk" <- function(e1, e2) {
  if (missing(e2)) {
    cli::cli_abort(c(
      "Cannot use {.code +} with a single argument",
      "i" = "Did you accidentally put {.code +} on a new line?"
    ))
  }

  # Get the name of what was passed in as e2, and pass along so that it
  # can be displayed in error messages
  e2name <- deparse(substitute(e2))

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

#' Add a `pk_data_settings` object.
#'
#' @param object The `pk_data_settings` object to be added.
#' @param pk_obj The `pk` object to which the `pk_data_settings` object will be added.
#' @param objectname The name of the `pk_data_settings` object.
#'
#' @return The `pk` object, modified by the `pk_data_settings` object.
#' @author Caroline Ring
#' @export
pk_add.pk_data_settings <- function(object, pk_obj, objectname){

  #New data_settings will *replace* existing ones
  if(!is.null(pk_obj$data_settings)){
    message(paste0(objectname,
                   ": data_settings already present; new data_settings will replace the existing one")
    )
  }

  pk_obj$data_settings <- object
  if(pk_obj$status > (status_preprocess - 1)){
    #with new data pre-processing settings, everything will change starting from
    #data pre-processing
    pk_obj$status <- status_preprocess - 1
    message(paste0(objectname,
                   ": Modifying data_settings resets status to level ",
                   status_preprocess - 1,
                   "; all later stages of the workflow will need to be re-done")
    )
  }
  return(pk_obj)
}

#' Add a `pk_optimx_settings` object.
#'
#' @param object The `pk_optimx_settings` object to be added.
#' @param pk_obj The `pk` object to which the `pk_optimx_settings` object will be added.
#' @param objectname The name of the `pk_optimx_settings` object.
#'
#' @return The `pk` object, modified by adding the settings.
#' @author Caroline Ring
#' @export
pk_add.pk_optimx_settings <- function(object, pk_obj, objectname){

  #New optimx_settings will *replace* existing ones
  if(!is.null(pk_obj$optimx_settings)){
    message(paste0(objectname,
                   ": optimx_settings already present; new optimx_settings will replace the existing one")
    )
  }

  pk_obj$optimx_settings <- object
  if(pk_obj$status > (status_prefit - 1)){
    #with new optimizer settings, data pre=processing and model pre-fitting
    #should not change, but model fitting will change
    message(paste0(objectname,
                   ": Modifying optimx_settings resets status to level ",
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
  #New stat_models will *replace* existing ones by the same name
  for(this_model in names(object)){
    if(!is.null(pk_obj$stat_model[[this_model]])){
      message(paste0(objectname,
                     ": stat_model for",
                     this_model,
                     "already present; new stat_model will replace the existing one")
      )
    }
    pk_obj$stat_model[[this_model]] <- object[[this_model]]
  }
  if(pk_obj$status > (status_prefit - 1)){
    message(paste0(objectname,
                   ": Modifying stat_model resets status to level",
                   status_prefit - 1,
                   "; all later stages of the workflow will need to be re-done")
    )
    pk_obj$status <- (status_prefit - 1) #data pre-processing won't change with addition of new model, but model pre-fit and fit will change
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


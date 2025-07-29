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
#' @export
#' @import cli
#' @author Caroline Ring
"+.pk" <- function(e1, e2) {
  # throw error if only one argument
  if (missing(e2)) {
    cli::cli_abort(c(
      "Cannot use {.code +} with a single argument",
      "i" = "Did you accidentally put {.code +} on a new line?"
    ))
  }

  # Get the name of what was passed in as e2, and pass along so that it
  # can be displayed in error messages
  e2name <- deparse(substitute(e2))

  # make sue
  if (is.pk(e1)) add_pk(e1, e2, e2name)
  else if (is.pkproto(e1)) {
    cli::cli_abort(c(
      "Cannot add {.cls pkproto} objects together",
      "i" = "Did you forget to add this object to a {.cls pk} object?"
    ))
  }
}

#' Add various `pkproto` objects to a `pk` object
#'
#' @param pk_obj The `pk` object
#' @param pkproto_obj The `pkproto` object to be added
#' @param objectname The name of the `pkproto` object to be added
#'
#' @return The `pk` object modified by the addition.
#' @export
add_pk <- function(pk_obj, pkproto_obj, objectname) {
  if (is.null(pkproto_obj)) return(pk_obj)

  p <- pk_add(pkproto_obj, pk_obj, objectname)
  p
}

#' Add a `pkproto` object to a `pk` object
#'
#' This is the S3 generic method.
#'
#' @param pkproto_obj The `pkproto` object to be added
#' @param pk_obj The `pk` object to which the `pkproto` object is to be added
#' @param objectname The object name
#' @return The `pk` object modified by the addition.
#'
#' @export
#' @seealso [pk_add.pk_scales()] for the method for adding `pk_scales` objects
#'   (from [scale_conc()] and [scale_time()]); [pk_add.pk_settings_preprocess()]
#'   for the method for adding `pk_settings_preprocess` objects (from
#'   [settings_preprocess()]); [pk_add.pk_nca_group()] for the method
#'   for adding `pk_nca_group` objects (from [settings_data_info()]);
#'   [pk_add.pk_settings_optimx()] for the method for adding
#'   `pk_settings_optimx` objects (from [settings_optimx()]);
#'   [pk_add.pk_stat_model()] for the method for adding `pk_stat_model` objects
#'   (from `stat_model()`)
pk_add <- function(pkproto_obj, pk_obj, objectname) {
  UseMethod("pk_add")
}

#' Add pkproto object default method
#' @param pkproto_obj The `pkproto` object to be added
#' @param pk_obj The `pk` object to which the `pkproto` object is to be added
#' @param objectname The object name
#' @return The `pk` object modified by the addition.
#' @export
pk_add.default <- function(pkproto_obj, pk_obj, objectname) {
  stop("No 'pk_add' method exists for object of class",
       pkproto_obj)
}

#' Subtract a `pkproto` object from a `pk` object
#'
#' This is the S3 generic method.
#'
#' @param pkproto_obj The `pkproto` object to be subtracted
#' @param pk_obj The `pk` object to which the `pkproto` object is to be subtracted
#' @param objectname The object name
#' @return The `pk` object modified by the subtraction.
#'
#' @export
#' @seealso
#'   [pk_subtract.pk_stat_model()] for the method for subtracting `pk_stat_model` objects
#'   (from `stat_model()`)
pk_subtract <- function(pkproto_obj, pk_obj, objectname) {
  UseMethod("pk_subtract")
}

#' Subtract pkproto object default method
#' @param pkproto_obj The `pkproto` object to be subtracted
#' @param pk_obj The `pk` object to which the `pkproto` object is to be subtracted
#' @param objectname The object name
#' @return The `pk` object modified by the addition.
#' @export
pk_subtract.default <- function(pkproto_obj, pk_obj, objectname) {
  stop("No 'pk_subtract' method exists for object of class",
       pkproto_obj)
}

#' Add a `pk_scales` object to a `pk` object.
#'
#' @param pkproto_obj The `pk_scales` object to be added.
#' @param pk_obj The `pk` object to which the `pk_scales` object will be added.
#' @param objectname The name of the `pk_scales` object.
#'
#' @return The `pk` object, modified by the `pk_scales` object.
#' @author Caroline Ring
#' @export
pk_add.pk_scales <- function(pkproto_obj, pk_obj, objectname) {

  pk_obj$scales[[pkproto_obj$name]] <- pkproto_obj$value

  if (pk_obj$status > (status_preprocess - 1L)) {
    # with new scaling, everything will change starting from data pre-processing
    # with new data pre-processing settings, everything will change starting from
    # data pre-processing
    pk_obj$status <- status_preprocess - 1L
    message(objectname,
            ": Modifying data scaling resets status to level ",
            status_preprocess - 1L,
            "; all later stages of the workflow will need to be re-done"
    )
    pk_obj$status <- status_preprocess - 1L
  }

  return(pk_obj)
}

#' Add a `pk_settings_preprocess` object.
#'
#' @param pkproto_obj The `pk_settings_preprocess` object to be added.
#' @param pk_obj The `pk` object to which the `pk_settings_preprocess` object will be added.
#' @param objectname The name of the `pk_settings_preprocess` object.
#'
#' @return The `pk` object, modified by the `pk_settings_preprocess` object.
#' @author Caroline Ring
#' @export
pk_add.pk_settings_preprocess <- function(pkproto_obj, pk_obj, objectname) {

  # New settings_preprocess will *replace* existing ones
  if (!is.null(pk_obj$pk_settings$preprocess)) {
    message(objectname,
            ": settings_preprocess already present; ",
            "new settings_preprocess will replace the existing one"
    )
  }

  pk_obj$pk_settings$preprocess <- pkproto_obj

  if (pk_obj$status > (status_preprocess - 1L)) {
    # with new data pre-processing settings, everything will change starting from
    # data pre-processing
    pk_obj$status <- status_preprocess - 1L
    message(objectname,
            ": Modifying settings_preprocess resets status to level ",
            status_preprocess - 1L,
            "; all later stages of the workflow will need to be re-done"
    )
  }
  return(pk_obj)
}

#' Add a `pk_nca_group` object.
#'
#' @param pkproto_obj The `pk_nca_group` object to be added.
#' @param pk_obj The `pk` object to which the `pk_nca_group` object will be added.
#' @param objectname The name of the `pk_nca_group` object.
#'
#' @return The `pk` object, modified by the `pk_nca_group` object.
#' @author Caroline Ring, Gilberto Padilla Mercado
#' @export
pk_add.pk_nca_group <- function(pkproto_obj, pk_obj, objectname) {

  # New settings_data_info will *replace* existing ones
  if (!is.null(pk_obj$pk_groups$nca_group)) {
    message(objectname,
            ": nca_group already present; ",
            "new nca_group will replace the existing one"
    )
  }

  pk_obj$pk_groups$nca_group <- pkproto_obj

  if (pk_obj$status > (status_preprocess - 1L)) {
    # with new data info settings, everything will change starting from
    # data info
    pk_obj$status <- status_data_info - 1L
    message(objectname,
            ": Modifying nca_group resets status to level ",
            status_data_info - 1L,
            "; all later stages of the workflow will need to be re-done"
    )
  }
  return(pk_obj)
}

#' Add a `pk_settings_optimx` object.
#'
#' @param pkproto_obj The `pk_settings_optimx` object to be added.
#' @param pk_obj The `pk` object to which the `pk_settings_optimx` object will be added.
#' @param objectname The name of the `pk_settings_optimx` object.
#'
#' @return The `pk` object, modified by adding the settings.
#' @author Caroline Ring
#' @export
pk_add.pk_settings_optimx <- function(pkproto_obj, pk_obj, objectname) {

  # New settings_optimx will *replace* existing ones
  if (!is.null(pk_obj$pk_settings$optimx)) {
    message(objectname,
            ": settings_optimx already present; ",
            "new settings_optimx will replace the existing one"
    )
  }

  pk_obj$pk_settings$optimx <- pkproto_obj

  if (pk_obj$status > (status_prefit - 1L)) {
    # with new optimizer settings, data pre=processing and model pre-fitting
    # should not change, but model fitting will change
    message(objectname,
            ": Modifying settings_optimx resets status to level ",
            (status_prefit - 1L),
            "; all later stages of the workflow will need to be re-done"
    )
    pk_obj$status <- 3L
  }
  return(pk_obj)
}


#' Add a `pk_stat_model` object.
#'
#' @param pkproto_obj The `pk_stat_model` object to be added.
#' @param pk_obj The `pk` object to which the `pk_stat_model` object will be added.
#' @param objectname The name of the `pk_stat_model` object.
#'
#' @return The `pk` object, modified by adding the `stat_model`.
#' @author Caroline Ring
#' @export
pk_add.pk_stat_model <- function(pkproto_obj, pk_obj, objectname) {
  # New stat_model will *replace* existing stat_model
  if (!is.null(pk_obj$stat_model)) {
    old_models <- names(pk_obj$stat_model)
    new_models <- names(pkproto_obj)
    message(objectname,
            ": stat_model already present; ",
            "new stat_model will replace the existing one.",
            "\n",
            "Old models: ", toString(old_models),
            "\n",
            "New models: ", toString(new_models)
    )
  }
  pk_obj$stat_model <- pkproto_obj

  # data pre-processing won't change with addition of new model, but model
  # pre-fit and fit will change
  if (pk_obj$status > (status_prefit - 1L)) {
    message(objectname,
            ": Modifying stat_model resets status to level ",
            status_prefit - 1L,
            "; all later stages of the workflow will need to be re-done"
    )
    pk_obj$status <- (status_prefit - 1L)
  }

  return(pk_obj)
}

#' Add a `pk_stat_error_model` object.
#'
#' @param pkproto_obj The `pk_stat_error_model` object to be added.
#' @param pk_obj The `pk` object to which the `pk_stat_error_model` object will be added.
#' @param objectname The name of the `pk_stat_error_model` object.
#'
#' @return The `pk` object, modified by adding the `stat_error_model`.
#' @author Caroline Ring
#' @export
pk_add.pk_stat_error_model <- function(pkproto_obj, pk_obj, objectname) {

  if (!is.null(pk_obj$stat_error_model)) {
    message(objectname,
            ": stat_error_model already present; ",
            "new stat_error_model will replace the existing one"
    )
  }

  pk_obj$pk_groups$error_group <- pkproto_obj

  # data pre-processing won't change with addition of a new error model, but
  # model pre-fit and fit will change
  if (pk_obj$status > (status_prefit - 1L)) {
    message(objectname,
            ": Modifying stat_error_model resets status to level ",
            status_prefit - 1L,
            "; all later stages of the workflow will need to be re-done"
            )
    pk_obj$status <- status_prefit - 1L
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
#' @param pkproto_obj The `uneval` (mapping) object to be added.
#' @param pk_obj The [pk()] object to which the `uneval` object will be added.
#' @param objectname The name of the `uneval` object.
#'
#' @return The [pk()] object, modified by adding the new mapping.
#' @author Caroline Ring
#' @export
pk_add.uneval <- function(pkproto_obj, pk_obj, objectname) {
  if (!is.null(pk_obj$mapping)) {
    message(objectname,
            ": mapping already present; new mapping will replace the existing one"
    )
  }

  pk_obj$mapping <- pkproto_obj

  # data pre-processing and everything downstream will change
  if (pk_obj$status > (status_preprocess - 1L)) {
    message(objectname,
            ": Modifying mapping resets status to level ",
            status_preprocess - 1L,
            "; all later stages of the workflow will need to be re-done"
            )
    pk_obj$status <- status_preprocess - 1L
  }

  return(pk_obj)
}

#' Add facet_data()
#'
#' @param pkproto_obj The `pk_facet_data` object to be added.
#' @param pk_obj The [pk()] object to which the `pk_facet_data` object will be added.
#' @param objectname The name of the `pk_facet_data` object.
#'
#' @return The [pk()] object, with the added `pk_facet_data` in the `groups` sub-list.
#' @author Caroline Ring
#' @export
pk_add.pk_facet_data <- function(pkproto_obj, pk_obj, objectname) {

  # data pre-processing and everything downstream will change
  if (pk_obj$status > (status_preprocess - 1L)) {
    message(objectname,
            ": Modifying faceting resets status to level ",
            status_preprocess - 1L,
            "; all later stages of the workflow will need to be re-done"
    )
    pk_obj$status <- status_preprocess - 1L
  }

  # this will become the data_group
  # this takes the harmonized variables in the pk object
    pk_obj$pk_groups$data_group <- pkproto_obj
    return(pk_obj)
}


#' Add loq_group
#'
#' @param pkproto_obj The `pk_loq_group` object to be added.
#' @param pk_obj The [pk()] object to which the `pk_loq_group` object will be added.
#' @param objectname The name of the `pk_loq_group` object.
#'
#' @return The [pk()] object, with the added `pk_loq_data` in the `groups` sub-list.
#' @author Gilberto Padilla Mercado
#' @export
pk_add.pk_loq_group <- function(pkproto_obj, pk_obj, objectname) {
  if (pk_obj$status > (status_preprocess - 1L)) {
    message(objectname,
            ": Modifying mapping resets status to level ",
            status_preprocess - 1L,
            "; all later stages of the workflow will need to be re-done"
    )
    pk_obj$status <- status_preprocess - 1L
  }

  # add loq to pk object 'groups' list
  pk_obj$pk_groups$loq_group <- pkproto_obj
  return(pk_obj)
}

#' Add sd_group
#'
#' @param pkproto_obj The `pk_sd_group` object to be added.
#' @param pk_obj The [pk()] object to which the `pk_sd_group` object will be added.
#' @param objectname The name of the `pk_sd_group` object.
#'
#' @return The [pk()] object, with the added `pk_sd_group` in the `groups` sub-list.
#' @author Gilberto Padilla Mercado
#' @export
pk_add.pk_sd_group <- function(pkproto_obj, pk_obj, objectname) {
  if (pk_obj$status > (status_preprocess - 1L)) {
    message(objectname,
            ": Modifying mapping resets status to level ",
            status_preprocess - 1L,
            "; all later stages of the workflow will need to be re-done"
    )
    pk_obj$status <- status_preprocess - 1L
  }

  # add loq to pk object 'groups' list
  pk_obj$pk_groups$sd_group <- pkproto_obj
  return(pk_obj)
}


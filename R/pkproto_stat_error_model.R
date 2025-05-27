#' Error model
#'
#' Define an error model.
#'
#' `stat_error_model` defines groupings for a fixed-effects error model. For each
#' model in `stat_model`, a single set of model parameters will be fit to `data`.
#' In order to do the fitting, the residual errors (observed concentrations -
#' model-predicted concentrations) are assumed to obey a zero-mean normal
#' distribution. However, in this package, the residuals are not all required to
#' obey the *same* zero-mean normal distribution. Different groups of residuals
#' may obey zero-mean normal distributions with different variances.
#' `stat_error_model` defines these groups as unique combinations of the
#' variables given in argument `error_group`. For example, the default value
#' `vars(Chemical, Species, Reference, Media)` means that for each group of
#' observations in `data` with a unique combination of `Chemical`, `Species`,
#' `Reference`, and `Media`, there is a separate residual error variance. For
#' example, if there happened to be three such unique combinations, there would
#' be three error variances.
#'
#' If you want all residuals to obey the same zero-mean normal distribution
#' (i.e., for there to be only one residual error variance), then you should
#' provide an `error_group` that puts all the data in the same group. For
#' example, since all data in `data` should already be for a single `Chemical`
#' and `Species`, you could provide `error_group = vars(Chemical, Species)` to
#' put all the data in the same group.
#'
#' Note that, since all data in `data` should already be for a single `Chemical`
#' and `Species`, you could leave out `Chemical` and `Species` from `error_group`
#' and still get the same result. However, we recommend explicitly including
#' `Chemical` and `Species`. Tncluding them will make your code more explicit and
#' transparent, and it does no harm. In addition, [invivopkfit] may be extended
#' in the future to allow input of data with multiple chemicals or species;
#' explicitly including `Chemical` and `Species` in your `error_group` will
#' future-proof your code in that sense.
#'
#' The error variance(s) are hyperparameters that will be estimated from the data
#' along with the model parameters. That means there needs to be enough data to
#' fit the model parameters plus the error variances. For example, if you are
#' fitting a 1-compartment model to oral and IV data measured in plasma, and
#' using an error model with three separate error-variance groups (e.g. three
#' different References), then you are trying to fit 4 model parameters (kelim,
#' Vdist, Fgutabs, kgutabs) plus 3 error variances, for a total of 7 parameters.
#' That means you need to have at least 8 data points. (When you call [prefit()],
#' this checking is done automatically. But it is useful to be aware of this, in
#' case you are trying to figure out why your fit was aborted due to insufficient data availability.)
#'
#' @param error_group Defined using [dplyr::vars()]: A set of harmonized
#'  variables whose unique combinations define a group with its own error
#'  variance. These variables refer to the `data` element of the `pk` object to
#'  which `stat_error_model()` is added. The variables should not be quoted.
#'  Default is `vars(Chemical, Species, Reference)`.
#' @param ... Additional arguments. Not currently used.
#' @return An object of class `pk_stat_error_model`: A named list of all the
#'  arguments to `stat_error_model`.
#' @author Caroline Ring
#' @export
stat_error_model <- function(error_group = vars(Chemical, Species, Reference),
                             ...) {
  # get arguments and values as a list
  argg <- c(as.list(environment()), list(...))
  this_error_model <- argg
  class(this_error_model) <- c(class(this_error_model),
                               "pkproto",
                               "pk_stat_error_model")
  return(this_error_model)
}

#' Setup data method generic
#'
#' For `pk` objects, see [setup_data.pk()].
#'
#' @param x An object.
setup_data <- function(x){
  UseMethod("setup_data", x)
}

setup_data.default <- function(x){
  stop(paste("No 'setup_data' method exists for object of class",
       class(x)))
}

#' Non-compartmental data method generic
#'
#' For `pk` objects, see [nca.pk()]
#'
#' @param x An object.
nca <- function(x){
  UseMethod("nca", x)
}

nca.default <- function(x){
  stop(paste("No 'nca' method exists for object of class",
             class(x)))
}

plot <- function(x){
  UseMethod("plot", x)
}

plot_data <- function(x){
  UseMethod("plot_data", x)
}

fit <- function(x){
  UseMethod("fit", x)
}

plot_fit <- function(x){
  UseMethod("plot_fit", x)
}

setup_analysis <- function(x){
  UseMethod("setup_analysis", x)
}

analyze <- function(x){
  UseMethod("analyze", x)
}

plot_analysis <- function(x){
  UseMethod("plot_analysis", x)
}

get_summary <- function(x){
  UseMethod("get_summary", x)
}

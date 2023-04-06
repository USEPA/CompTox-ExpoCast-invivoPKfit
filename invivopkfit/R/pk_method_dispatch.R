setup_data <- function(x){
  UseMethod("setup_data", x)
}

nca <- function(x){
  UseMethod("nca", x)
}

plot <- function(x){
  UseMethod("plot", x)
}

plot_data <- function(x){
  UseMethod("plot_data", x)
}

setup_fit <- function(x){
  UseMethod("setup_fit", x)
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

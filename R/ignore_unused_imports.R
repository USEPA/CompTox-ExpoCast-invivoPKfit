#' Ignore unused imports
#'
#' Placeholder function to appease R CMD CHECK
#'
#' This function does nothing and should be ignored by the user.
#'
#' Why it exists: Whenever possible, `invivopkfit` code calls functions from other packages
#' using the syntax `package::function()`, which means that the whole package
#' does not have to be loaded, nor does the function itself have to be loaded
#' until it is used. The relevant packages are listed in the `invivofit` package
#' `DESCRIPTION` file under `Imports`, because they must be installed to use
#' `invivopkfit`. But because the packages are not actually loaded, this creates
#' a (spurious) NOTE from `R CMD CHECK` about a declared Import not being used.
#' This function is a workaround to suppress that NOTE. It does nothing except contain
#' namespace-qualified references (not calls) to objects in the relevant packages.
#'
#' @author Caroline Ring
#' @references https://r-pkgs.org/dependencies-in-practice.html#how-to-not-use-a-package-in-imports
#'
ignore_unused_imports <- function() {
  pracma::cumtrapz

}

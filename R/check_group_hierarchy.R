#' Checking data, error, and summary group hierarchical structure
#'
#' @param obj An object created by [pk()].
#'
#' @returns Prints a tree of the summary hierarchies and an error when
#'   the hierarchical structure expectation is not met.
#'
check_group_hierarchy <- function(obj) {
  data_group <- get_data_group.pk(obj, as_character = TRUE)
  error_group <- get_error_group.pk(obj, as_character = TRUE)
  summary_group <- get_nca_group.pk(obj, as_character = TRUE)

  ehg_check <- all(data_group %in% error_group)
  shg_check0 <- all(error_group %in% summary_group)
  shg_check1 <- all(data_group %in% summary_group)

  dh_group <- paste(data_group, collapse = "--")
  eh_group <- paste(base::setdiff(error_group, data_group), collapse = "--")
  sh_group <- paste(base::setdiff(summary_group, error_group), collapse = "--")

  tree_df <- data.frame(
    group = c(dh_group, eh_group, sh_group),
    subgroup = I(list(eh_group, sh_group, character(0L))),
    group_name = c(
      paste0("Data Group: [", dh_group, "]"),
      paste0(" + Error Group: [", eh_group, "]"),
      paste0(" + Summary Group: [", sh_group, "]")
    ),
    stringsAsFactors = FALSE
  )

  if (!all(ehg_check, shg_check0)) {
    cli::cli_abort(
      c(
        "Hierarchy expectation in model is not met.",
        "!" = paste(
          "Please check that all data group variables are a subset of error group variables",
          "and all error group variables are a subset of summary group variables."
        )
      )
    )
  }

  cli::cli_inform(c("v" = "Current Modeling Structure:"))
  print(cli::tree(data = tree_df))
  cli::cli_inform("")

}

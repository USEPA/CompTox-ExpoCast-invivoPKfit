pk_methods_list <- dir(path = "R", pattern = "pk_methods")

pk_methods_names <- gsub(x=pk_methods_list,
                         pattern = "pk_methods_",
                         replacement = "",
                         fixed = TRUE)
pk_methods_names <- gsub(x=pk_methods_names,
                         pattern = ".R",
                         replacement = "",
                         fixed = TRUE)

#these behave differently (they modify in place)
#and I have created these manually.
pk_methods_names <- base::setdiff(pk_methods_names,
                            c("preprocess_data",
                              "data_info",
                              "prefit",
                              "fit",
                              "add_pkproto"))

for (i in seq_along(pk_methods_names)){
  this_name <- pk_methods_names[i]

  #produce text of the pk_faceted script

  this_text <- c(paste0("#'",
                        c(this_name,
                          "",
                          paste(this_name,
                                "for pk_faceted objects"),
                          "",
                          "@param obj An object of class `pk_faceted`",
                          paste("@return",
                                "A [tibble::tibble()] grouped by the faceting variables",
                                "with a new variable `",
                                this_name,
                                "` containing the results of applying",
                                paste0("[", this_name, ".pk()]"),
                                "to each item in the list column of [pk()] objects, `obj$pk_object`"),
                          "@export",
                          "@author Caroline Ring",
                          "@family methods for pk_faceted objects"
                          )
                        ),
                 paste0(this_name, ".pk_faceted <- function(obj, ...){"),
                 paste0("obj %>% dplyr::mutate(",
                        this_name, " =",
                        "purrr::map(pk_object,",
                        this_name, "))"),
                 "}"
                 )

  writeLines(this_text,
             paste0("R/pk_faceted_methods_", this_name, ".R"))

}

rename_vars <- function(data, ...){
  dots <- enquos(...)
  as.data.frame(sapply(dots,
         function(q) eval_tidy(q, data),
         simplify = FALSE,
         USE.NAMES = TRUE))
}

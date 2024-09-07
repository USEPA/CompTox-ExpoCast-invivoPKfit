stat_nca <- function(summary_group = dplyr::vars(Chemical,
                                             Species,
                                             Reference,
                                             Route,
                                             Media,
                                             Dose),
                     ...){
  #get arguments and values
  argg <- c(as.list(environment()), list(...))
  this_stat_nca <- argg
  #set class
  class(this_stat_nca) <- c(class(this_stat_nca),
                            "pkproto",
                            "pk_stat_nca")

  return(this_stat_nca)
}

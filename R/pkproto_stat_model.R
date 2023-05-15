#' Add PK model(s) to be fitted
#'
#' @param model The name(s) of models to be fitted. These should be the names of objects of class `pk_model`. Built-in options are `"flat"`,
#'   `"1comp"`, and `"2comp"`. You may add your own model by using [pk_model()].

stat_model <- function(model = c("flat", "1comp", "2comp"),
                      ...){
  this_stat_model <- list()
  #pull pk_model objects by each name
for(this_model in model){
  this_stat_model[[this_model]] <- list()
  #check whether an object exists by this name
  if(exists(this_model)){
    this_model_obj <- get(this_model)
    #Check whether this is an object of class `pk_model`
    if(inherits(this_model_obj, "pk_model")){
      #if so, add it to the list
      this_stat_model[[this_model]] <- this_model_obj
    }
  }
}

  #set class
  class(this_stat_model) <- c(class(this_stat_model), "pkproto", "pk_stat_model")

  return(this_stat_model)
}


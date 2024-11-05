paste2 <- function(...,
                   sep = " ",
                   collapse = NULL,
                   na.rm = TRUE,
                   trimws = TRUE,
                   blank.rm = TRUE){
  arglist <- list(...)
  #convert the arguments to character
  arglist <- lapply(arglist,
                    as.character)
  #trim whitespace on the arguments
  if(trimws %in% TRUE){
    arglist <- lapply(arglist, trimws)
  }

  #set blanks to NA if so specified
  if(blank.rm %in% TRUE){
    arglist <- lapply(arglist,
                      function(x){
                        x[length(x) == 0] <- NA_character_
                        x[!nzchar(x) %in% TRUE] <- NA_character_
                        x
                      })
  }
  #cbind the arguments
  argc <- do.call(cbind, arglist)
  #do sep-phase
  out <- apply(argc,
        1,
        function(x) paste(x[!is.na(x)],
                          collapse = sep))

  #do collapse phase
  if(!is.null(collapse)){
  out <- paste2(out,
                        sep = collapse,
                        collapse = NULL,
                        na.rm = na.rm,
                        trimws = trimws,
                        blank.rm = blank.rm)

  }

  return(out)
}



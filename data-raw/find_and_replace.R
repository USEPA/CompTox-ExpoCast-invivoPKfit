find_and_replace <- function(path = "R",
                             pattern= NULL,
                             find = NULL,
                             replace = NULL){

  flist <- dir(path = path,
               pattern= pattern,
               full.names = TRUE)

  junk <- sapply(flist,
         function(this_file){
           tmp <- readLines(con = this_file)

           tmp <- gsub(pattern = find,
                       replacement = replace,
                       x = tmp)

           writeLines(tmp,
                      con = this_file)

         })

  return(0)
}







#' Combined mean
#'
#' Given mean, standard deviation, and N for some set of groups, calculate the
#' combined mean. Note that the groups may not overlap.
#'
#' @param group_mean Numeric vector: Observed sample means for summary data, or
#'   observed values for non-summary data. Censored observations should *not* be
#'   NA; they should be substituted with some value at or below the
#'   corresponding LOQ (e.g. LOQ or LOQ/2). Even if `log %in% TRUE`, these
#'   should *not* be log-transformed.
#' @param group_n Numeric vector: Observed sample number of subjects for summary
#'   data. For non-summary data (individual-subject observations), `group_n`
#'   should be set to 1.
#' @param na.rm Logical. If TRUE (default), then any groups where mean *or*
#'   N were NA will be dropped. If FALSE, they will be retained (and the result
#'   will be NA).
#' @return Numeric: the mean of the combined population (i.e. if
#'   all the groups were concatenated into one large group).
#' @author Caroline Ring
combined_mean <- function(group_mean,
                          group_n,
                          na.rm = TRUE,
                          log = FALSE){

  x_len <- c("group_mean" = length(group_mean),
             "group_n" = length(group_n))

  if(any(x_len %in% 0)){
    stop(paste0("invivopkfit::combined_mean(): ",
                "the following arguments have zero length: ",
                paste(names(x_len)[x_len %in% 0],
                      collapse = ", ")
    ))
  }

  max_len <- max(x_len)
  which_max_len <- which.max(x_len)

  bad_len <- (x_len < max_len) & (x_len != 1)


  if(any(bad_len)){
    warning(paste("invivopkfit::combined_sd():",
                  "the following inputs do not have matching lengths: ",
                  paste(paste0(names(x_len)[bad_len],
                               " length = ",
                               x_len[bad_len]),
                        collapse = "\n"
                  ),
                  "\n They will be repeated to match the length of the longest input,",
                  names(x_len)[which_max_len],
                  " length = ",
                  max_len,
                  "."
    ))
  }

  #repeat to match longest
  for (i in seq_along(x_len)){
    assign(names(x_len)[i],
           rep( #repeat the current value of each item to match the length
             get(names(x_len)[i]), #get the current value of each item
             length.out = max_len)
    )
  }

  #remove NAs if so specified
  if(na.rm %in% TRUE){
    which_na <- is.na(group_mean) | is.na(group_sd) | is.na(group_n)
    group_mean <- group_mean[!which_na]
    group_n <- group_n[!which_na]
  }

  grand_mean <- sum(group_n*group_mean)/sum(group_n)

  return(grand_mean)
}

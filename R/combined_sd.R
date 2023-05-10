#' Combined standard deviation
#'
#' Given mean, standard deviation, and N for some set of groups, calculate the
#' combined standard deviation. Note that the groups may not overlap.
#'
#' @param group_mean Numeric vector: Observed sample means for summary data, or
#'   observed values for non-summary data. Censored observations should *not* be
#'   NA; they should be substituted with some value at or below the
#'   corresponding LOQ (e.g. LOQ or LOQ/2). Even if `log %in% TRUE`, these
#'   should *not* be log-transformed.
#' @param group_sd Numeric vector: Observed sample SDs for summary data. For
#'   non-summary data (individual-subject observations), the corresponding
#'   element of `group_sd` should be set to 0. Even if `log %in% TRUE`, these
#'   should *not* be log-transformed.
#' @param group_n Numeric vector: Observed sample number of subjects for summary
#'   data. For non-summary data (individual-subject observations), `group_n`
#'   should be set to 1.
#' @param unbiased Logical. If TRUE (the default), then `group_sd` is assumed to
#'   be the unbiased estimator of population standard deviation (i.e. calculated
#'   using `n-1` in the denominator -- the way that `stats::sd()` calculates
#'   it), and the returned combined SD is also the unbiased estimator of the
#'   combined population SD. If FALSE, then `group_sd` is assumed to be the
#'   biased estimator (using `n` in the denominator), and the returned value is
#'   also the biased estimator of the combined population SD.
#' @param na.rm Logical. If TRUE (default), then any groups where mean, SD, *or*
#'   N were NA will be dropped. If FALSE, they will be retained (and the result
#'   will be NA).
#' @return Numeric: the standard deviation of the combined population (i.e. if
#'   all the groups were concatenated into one large group).
#' @author Caroline Ring
combined_sd <- function(group_mean,
                        group_sd,
                        group_n,
                        unbiased = TRUE,
                        na.rm = TRUE,
                        log = FALSE){

  x_len <- c("group_mean" = length(group_mean),
             "group_sd" = length(group_sd),
             "group_n" = length(group_n))

  if(any(x_len %in% 0)){
    stop(paste0("invivopkfit::combined_sd(): ",
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
    group_sd <- group_sd[!which_na]
    group_n <- group_n[!which_na]
  }

  grand_mean <- sum(group_n*group_mean)/sum(group_n)

  #if all N = 1, then just take regular standard deviation
  if(all(group_n %in% 1)){
    grand_sd <- sd(group_mean)

    grand_n <- sum(group_n)
    if(unbiased %in% FALSE){
      #convert unbiased SD to biased SD
      grand_sd <- grand_sd * sqrt((grand_n-1)/grand_n)
    }
  }else{ #if not all N = 1
    if(unbiased %in% TRUE){
      #convert unbiased group SDs to biased group SDs
      group_sd <- group_sd *
        sqrt((group_n[which_na]-1)/group_n)
    }

    grand_var <- (sum(group_n*group_sd^2) +
                    sum(group_n*(group_mean - grand_mean)^2))/
      (sum(group_n))

    grand_sd <- sqrt(grand_var) #biased


    if(unbiased %in% TRUE){
      grand_n <- sum(group_n)
      #convert biased grand SD to unbiased grand SD
      grand_sd <- grand_sd * sqrt(grand_n/(grand_n - 1))
    }
  }

  if(log %in% TRUE){
    #convert to log-scale combined SD
    grand_sd <- sqrt(log(1 + grand_sd^2 / grand_mean^2))
  }

  return(grand_sd)
}

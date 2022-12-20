#' Normal distribution density function for summary data
#'
#' Evaluates the normal distribution density function for summary data reported
#' as sample mean, sample SD, and sample N.
#'
#' `x_mean`, `x_sd`, and `X_N` must either be length 1, or the same size as the
#' others. If two or more are different lengths and not length 1, the function
#' will stop with an error.
#'
#' If `mu` and `sigma` are not the same
#'
#' @param mu Mean of the normal distribution to be evaluated (*not* the sample
#'   mean). May be a numeric scalar or vector.
#' @param sigma Standard deviation of the normal distribution to be evaluated
#'   (*not* the sample SD). May be a numeric scalar or vector.
#' @param x_mean Sample mean. May be a numeric scalar or vector.
#' @param x_sd Sample standard deviation. May be a numeric scalar or vector.
#' @param x_N Sample number of observations. May be a numeric scalar or vector.
#' @param log TRUE/FALSE: Whether to return the log of the density function.
#'   Default FALSE (to return the density function value on the natural scale).
#'
#' @return A numeric scalar or vector matching the length of the longest of
#'   `mu`, `sigma`, `x_mean`, `x_sd`, and `x_N`.
#' @export

dnorm_summary <- function(mu,
                             sigma,
                             x_mean,
                             x_sd,
                             x_N,
                             log = FALSE){

  x_len <- c("x_mean" = length(x_mean),
             "x_sd" = length(x_sd),
             "x_N" = length(x_N))


  if(sum(x_len %in% 1) == 2){
    #if two of them are length 1, repeat to match the longer
    if(x_len["x_mean"] %in% 1) x_mean <- rep(x_mean, length.out = max(x_len))
    if(x_len["x_sd"] %in% 1) x_sd <- rep(x_sd, length.out = max(x_len))
    if(x_len["x_N"] %in% 1) x_N <- rep(x_N, length.out = max(x_len))
  }else if(sum(x_len %in% 1) == 1 &
     length(unique(x_len)) == 2 ){
    #if one of them is length 1 and the others have the same length,
    #repeat to match the longest
    if(x_len["x_mean"] %in% 1) x_mean <- rep(x_mean, length.out = max(x_len))
    if(x_len["x_sd"] %in% 1) x_sd <- rep(x_sd, length.out = max(x_len))
    if(x_len["x_N"] %in% 1) x_N <- rep(x_N, length.out = max(x_len))
  }else if(sum(x_len %in% 1)==0 &
           length(unique(x_len)) > 1){
    #if none are length 1 and they have different lengths, throw error
    stop("invivopkfit::dnorm_summary(): x_mean, x_sd, and x_N must either be all the same length, or length 1.")
  }else if(sum(x_len %in% 1)==1 &
           length(unique(x_len)) > 2){
    #if one is length 1 and the other 2 have different lengths, throw error
    stop("invivopkfit::dnorm_summary(): x_mean, x_sd, and x_N must either be all the same length, or length 1.")
  }

  #Evaluate
  y_log <- x_N * log(1/(sigma*sqrt(2*pi))) +
    (-1/(2*sigma^2)) * (
      (x_N - 1) * x_sd^2 +
                         x_N * (x_mean - mu)^2
      )

 if(log == TRUE){
   return(y_log)
 }else{
   return(exp(y_log))
 }

}

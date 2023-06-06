sigma <- function(t, J, time, conc, conc_sd, n_subj) {
  sigma <- matrix(nrow = 2 * J, ncol = 2 * J)
  sigma[(1:(2 * J)), (1:(2 * J))] <- 0
  data <- cbind(conc, log(conc), time)
  #variance of concentrations for each unique time
  sx <- tapply(data[, 1], data[, 3], var)
  #variance of log concentrations for each unique time
  sy <- tapply(data[, 2], data[, 3], var)
  #loop over unique times:
  for (i in 1:J) {
    #get data for this unique time
    dintern <- subset(data, data[, 3] == t[i])
    #covariance of conc & log conc for this unique time
    #assigned to row i + J, where J is the number of unique time points,
    #and column i
    #lower triangular part of the matrix?
    sigma[i + J, i] <- cov(dintern[, 1], dintern[, 2])
    #and the upper triangular part of the matrix
    sigma[i, i + J] <- sigma[i + J, i]
  }
  #variances on the diagonal
  diag(sigma) <- c(sx, sy)
  sigma <- ifelse(is.na(sigma), 0, sigma)
  return(sigma)
}

pkparms <- function(time, conc, conc_sd, n_subj, dose, n.tail) {

  k <- length(unique(time)) - n.tail #keep 1:k points
  #trapezoidal rule weights
  w <- PK:::.weight(unique(time))

  t <- unique(time)
  J <- length(t)
  d <- dose

  M <- sigma(t = t, J = J, time = time, conc = conc)

  x <- conc
  x[x == 0] <- 1e-256

  df <- data.frame(time = time,
                   conc = conc,
                   conc_sd = conc_sd,
                   n_subj = n_subj)
  df_log <- as.data.frame(convert_summary_to_log10(sample_mean = df$conc,
                                     sample_sd = df$conc_sd))
  df_log$n_subj <- df$n_subj

  #average conc for each unique time
  xq <- df %>% dplyr::group_by(time) %>%
    dplyr::summarise(grand_mean = combined_mean(group_mean = conc,
                                                group_n = n_subj)) %>%
    dplyr::pull(grand_mean)

  #variance for each unique time
  sx <- df %>% dplyr::group_by(time) %>%
    dplyr::summarise(grand_sd = combined_sd(group_mean = conc,
                                            group_sd = conc_sd,
                                            group_n = n_subj)) %>%
    dplyr::pull(grand_sd)

  # xq <- as.vector(tapply(x,
  #                        time,
  #                        mean))

  y <- df_log$log10mean   #log scale concs
  #average log-scale conc for each unique time
  yq <- df_log %>% dplyr::group_by(time) %>%
    dplyr::summarise(grand_mean = combined_mean(group_mean = log10mean,
                                                group_n = n_subj)) %>%
    dplyr::pull(grand_mean)

  #log-scale variance for each unique time
  sy <- df_log %>% dplyr::group_by(time) %>%
    dplyr::summarise(grand_sd = combined_sd(group_mean = log10mean,
                                            group_sd = log10SD,
                                                group_n = n_subj)) %>%
    dplyr::pull(grand_sd)
  #yq <- tapply(y, time, mean, na.rm = TRUE) #avg log conc for each unique time
  #sy <- tapply(y, time, var, na.rm = TRUE) #variance of log conc for each unique time

  i <- c((k + 1):J) #index of points in the tail

  a <- c(rep(NA, k), t[i] - mean(t[i]))

  u <- c(rep(NA, k), a[i]/sum(a[i]^2))

  lambda <- sum(u[i] * yq[i]) * (-1)

  u <- u[!is.na(u)]
  i <- c(1:(J - 1))
  vec <- c(w[i], w[J] + 1/lambda, rep(0, k), xq[J] * u/lambda^2)
  est <- sum(w * xq) + xq[J]/lambda
  var <- (vec %*% M %*% vec)
  auc <- c(est, var)

  #AUC to tlast
  est <- sum(w * xq)
  vec <- c(w, rep(0, k + n.tail))
  var <- (vec %*% M %*% vec)
  auc_tlast <- c(est, var)

  #AUMC
  vec <- c(t[i] * w[i], t[J] * w[J] + t[J]/lambda + 1/lambda^2,
           rep(0, k), ((t[J] * lambda + 2) * xq[J] * u)/lambda^3)
  est <- sum(t * w * xq) + (xq[J]/lambda) * (t[J] + 1/lambda)
  var <- (vec %*% M %*% vec)
  aumc <- c(est, var)

  #MRT
  c1 <- (-aumc[1] - w[J] * lambda * aumc[1])/(lambda *
                                                auc[1]^2)
  c2 <- (lambda * t[J] + lambda^2 * t[J] * w[J] + 1)/(lambda^2 *
                                                        auc[1])
  c3 <- (2 * auc[1] - lambda * aumc[1] + lambda * t[J] *
           auc[1])/(lambda^3 * auc[1]^2)
  vec <- c(w[i] * (t[i]/auc[1] - aumc[1]/auc[1]^2), c1 +
             c2, rep(0, k), c3 * xq[J] * u)
  est <- aumc[1]/auc[1]
  var <- (vec %*% M %*% vec)
  mrt <- c(est, var)

  #Halflife
  hl <- c(log(2) * mrt[1], log(2)^2 * mrt[2])

  #Clearance
  vec <- c(-w[i] * d/auc[1]^2,
           (-d - d * lambda * w[J])/(lambda *
                                       auc[1]^2),
           rep(0, k),
           (-1) * (xq[J] * u * d)/(lambda^2 *
                                     auc[1]^2))
  est <- d/auc[1]
  var <- (vec %*% M %*% vec)
  cls <- c(est, var)

  #Vss
  A <- 2 * xq[J] + 2 * lambda * xq[J] * t[J] + 2 * lambda^2 *
    sum(t * w * xq)
  B1 <- auc[1] * lambda - 2 * xq[J] - 2 * lambda * xq[J] *
    t[J] - 2 * lambda * xq[J] * w[J] + auc[1] * lambda^2 *
    t[J]
  B2 <- auc[1] * lambda^3 * t[J] * w[J] - 2 * lambda^2 *
    sum(t[i] * w[i] * xq[i]) - 4 * lambda^2 * xq[J] *
    t[J] * w[J]
  B3 <- -2 * w[J] * lambda^3 * sum(t * w * xq)
  B <- B1 + B2 + B3
  C <- 2 * auc[1] * lambda - 2 * xq[J] - 2 * lambda * xq[J] *
    t[J] + auc[1] * lambda^2 * t[J] - 2 * lambda^2 *
    sum(t * w * xq)
  c1 <- (d * w[i] * t[i])/auc[1]^2 - (d * w[i] * A)/(auc[1]^3 *
                                                       lambda^2)
  c2 <- (d * B)/(auc[1]^3 * lambda^3)
  c3 <- (d * xq[J] * u * C)/(lambda^4 * auc[1]^3)
  vec <- c(c1, c2, rep(0, k), c3)
  est <- d * aumc[1]/auc[1]^2
  var <- (vec %*% M %*% vec)
  vss <- c(est, var)

  res <- rbind(auc, aumc, mrt, hl, cls, vss)
  rownames(res) <- c("AUC to infinity",  "AUMC to infinity",
                     "Mean residence time", "Half-life", "Clearance",
                     "Volume of Distribution")
  return(res)
}


# get.confint <- function(conc,
#                         time,
#                         k = k,
#                         d = d,
#                         alpha = alpha,
#                         nsample = 0) {
#   n <- mean(tapply(conc, time, length))
#   w <- PK:::.weight(unique(time))
#   t <- unique(time)
#   J <- length(t)
#   obsv.parms <- pkparms(time = time, conc = conc, d = d,
#                         k = k, w = w, t = t, J = J)
#   z <- qnorm(1 - alpha/2)
#   obsv.stderr <- sqrt(obsv.parms[, 2]/n)
#   asymp.lower <- obsv.parms[, 1] - obsv.stderr * z
#   asymp.upper <- obsv.parms[, 1] + obsv.stderr * z
#   asymp <- data.frame(est = obsv.parms[, 1],
#                       stderr = obsv.stderr,
#                       lower = asymp.lower,
#                       upper = asymp.upper,
#                       method = rep("z",
#                                    6),
#                       stringsAsFactors = TRUE)
#   res <- asymp
#   if (nsample > 0) {
#     boot.stat <- matrix(nrow = nsample, ncol = 6)
#     for (i in 1:nsample) {
#       boot.conc <- as.vector(unlist(tapply(conc,
#                                            time,
#                                            function(x){
#                                              if(length(x)==1){
#                                                x
#                                              }else{
#                                                sample(x, size = length(x))
#                                              }
#                                             },
#                                            replace = TRUE)
#                                     )
#                              )
#       boot.parms <- pkparms(conc = boot.conc, time = time,
#                             dose = dose, n.tail = n.tail)
#       boot.stat[i, ] <- (boot.parms[, 1] - obsv.parms[,
#                                                       1])/sqrt(boot.parms[, 2]/n)
#     }
#     t.lb <- apply(boot.stat, 2, quantile, probs = c(alpha/2),
#                   method = 5, na.rm = TRUE)
#     t.ub <- apply(boot.stat, 2, quantile, probs = c(1 -
#                                                       alpha/2), method = 5, na.rm = TRUE)
#     boott.lower <- obsv.parms[, 1] - t.ub * obsv.stderr
#     boott.upper <- obsv.parms[, 1] - t.lb * obsv.stderr
#     res <- data.frame(est = rep(asymp$est, each = 2),
#                       stderr = rep(asymp$stderr, each = 2),
#                       lower = as.vector(rbind(asymp.lower,
#                                               boott.lower)),
#                       upper = as.vector(rbind(asymp.upper,
#                                               boott.upper)),
#                       method = rep(c("z", "boott"),
#                                    6),
#                       stringsAsFactors = TRUE)
#   }
#   return(res)
# }
#
# nca.ssd <- function (conc, time, n.tail = 3, dose = 0, method = c("z", "boott"),
#           conf.level = 0.95, nsample = 1000, data)
# {
#
#   if (!missing(data)) {
#     cnames <- colnames(data)
#     if (!any(cnames == "conc")) {
#       stop("data does not contain a variable conc")
#     }
#     if (!any(cnames == "time")) {
#       stop("data does not contain a variable time")
#     }
#     conc <- data$conc
#     time <- data$time
#     if (any(cnames == "group")) {
#       group <- data$group
#     }
#     else {
#       group <- NULL
#     }
#   }
#   if (!is.vector(time) || !is.vector(conc)) {
#     stop("Argument time and/or conc invalid")
#   }
#   method <- match.arg(method, several.ok = TRUE)
#   if (!(any(method == "boott")))
#     nsample <- 0
#   if (any(method == "boott") & nsample == 0)
#     stop("boott called with zero resamples")
#   data <- data.frame(conc, time)
#   data <- data[order(data$time), ]
#   data <- na.omit(data)
#   data <- subset(data, data$conc >= 0)
#   conc <- data$conc
#   time <- data$time
#   k <- length(unique(time)) - n.tail
#   alpha <- 1 - conf.level
#   res <- get.confint(conc = conc, time = time, k = k, d = dose,
#                      alpha = alpha, nsample = nsample)
#   auc.obs <- PK::auc.ssd(conc = conc, time = time, method = levels(res$method),
#                      conf.level = conf.level, nsample = nsample)$CIs
#   if (any(method == "boott")) {
#     res <- rbind(auc.obs[2:1, -5], res)
#   }
#   else {
#     res <- rbind(auc.obs[1, -5], res)
#   }
#   r.names <- c("AUC to tlast", "AUC to infinity", "AUMC to infinity",
#                "Mean residence time", "non-compartmental half-life",
#                "Clearance", "Volume of distribution at steady state")
#   rownames(res) <- paste(conf.level * 100, "% CI for the ",
#                          rep(r.names, each = length(levels(res$method))), " using a ",
#                          sort(levels(res$method), decreasing = TRUE), " distribution",
#                          sep = "")
#   if (!any(method == "z")) {
#     res <- res[seq(2, 14, 2), ]
#   }
#   out <- NULL
#   out$CIs <- res
#   out$design <- "ssd"
#   out$est <- matrix(split(res, res$method)[[1]][, 1], ncol = 1)
#   rownames(out$est) <- r.names
#   colnames(out$est) <- "est"
#   out$conf.level <- conf.level
#   class(out) <- "PK"
#   out$conc <- conc
#   out$time <- time
#   out$group <- NULL
#   out$dose <- dose
#   return(out)
# }

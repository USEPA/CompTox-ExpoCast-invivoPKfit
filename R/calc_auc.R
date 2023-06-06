calc_auc <- function(time,
                     conc,
                     conc_sd,
                     n_subj,
                     ntail = 3){
data <- data.frame(time = time,
                   conc = conc,
                   conc_sd = conc_sd,
                   n_subj = n_subj)
#sort time from lowest to highest
data <- data[order(data$time), ]
#get unique time points
ut <- unique(data$time)
#get average concentration at each time
uc <- approx(x = data$time,
             y = data$conc,
             xout = ut,
             ties = "mean")$y
#get combined SD at each time
uc_sd <- data %>%
  dplyr::group_by(time) %>%
  dplyr::summarise(grand_sd = combined_sd(group_mean = conc,
                               group_sd = conc_sd,
                               group_n = n_subj)) %>%
  dplyr::pull(grand_sd)

if(!(min(ut)==0)){
  #if no zero time, add a zero time
  ut <- c(0, ut)
  #if oral, then extrapolate conc = 0 at time = 0
  if(route %in% "oral"){
    uc <- c(.Machine$double.eps, uc)
  }else{
    #if IV, then extrapolate back based on linear regression of first ntail points
    uc_first_log <- log(uc[1:ntail])
    ut_first <- ut[1:ntail]
    lm_first <- lm(uc_first_log ~ ut_first)
    uc <- c(predict(lm_first, newdata = 0),
            uc)
  }
}

#get AUC to tlast
auc_tlast <- pracma::trapz(x = ut,
              y = uc)

#propagation of uncertanity
#get trapezoidal rule weights
# tw <- diff(ut)/2

# #Taylor expansion
# del_auc_tlast <- c(tw, 0) + c(0, tw)
# var_auc_tlast <- del_auc_tlast^2 * uc_sd^2

#AUC to infinity: tail extrapolation
#slope of last n.tail points on a semilog scale
uc_last <- rev(uc)[1:ntail]
ut_last <- rev(ut)[1:ntail]
uc_last_log <- log(uc_last)
k <- -coef(lm(uc_last_log ~ ut_last))[2]
auc_inf <- auc_tlast + uc_last[1]/k

# #propagation of uncertanity
# #Taylor expansion
# del_auc_inf <- c(tw, 1/k) + c(0, tw)
# var_auc_inf <- del_auc_inf^2 * uc_sd^2

#AUMC
auc_cum <-  pracma::cumtrapz(x = ut,
                             y = uc)
auc_cum_last <- rev(auc_cum)[1:ntail]
y_last <- auc_inf - auc_cum_last
y_last_log <- log(y_last)
k_aumc <- -coef(lm(y_last_log ~ ut_last))[2]
aumc_inf <- pracma::trapz(x = ut,
                          y = auc_cum) + auc_cum_last[1]/k_aumc





}

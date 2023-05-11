calc_nca <- function(time,
                    dose,
                    conc,
                    detect,
                    n_subj,
                   subject_id,
                   n.tail = 3){
  #Calculate area under the concentration-time curve for NCA.
  #At present we do not know individual animal IDs or animal-group IDs for each point,
  #so we will have to just assume that each time point is a different animal.
  #substitute any nondetects with LOQ/2. This is not sophisticated but it is fast.

  conc <- ifelse(detect %in% FALSE,
                 conc*0.5,
                 conc)
  if(all(is.na(subject_id))){
    #assume every obs is a different subject
    subject_id <- rep(1L, length(conc))
  }


  ntab <- table(time, subject_id)
  m <- matrix(ntab, ncol = ncol(ntab))
  if(ncol(m)==1){
    design <- "complete"
  }else{
    if(all(m[!diag(nrow(m))] == 0)){
      #one measurement per subject per time
      design <- "ssd"
    }else{
      if(all(m==1)){
        design <- "complete"
      }else{
        design <- "batch"
      }
    }
  }

  pk_out <- tryCatch(PK::nca(conc = conc,
                   time = time,
                   n.tail = n.tail,
                   dose = unique(dose),
          design = design)$est[,1],
          error = function(err){
            tmp <- rep(NA_real_, 7)
            names(tmp) <- c("AUC to tlast",
                               "AUC to infinity",
                               "AUMC to infinity",
                               "Mean residence time",
                               "non-compartmental half-life",
                               "Clearance",
                               "Volume of distribution at steady state")
            tmp
          })

  #also compute tmax, Cmax
  peak <- unlist(get_peak(x = time,
                          y = conc))
  names(peak) <- c("tmax", "Cmax")
  pk_out <- c(pk_out, peak)

  outval <- data.frame("design" = design,
    "param_name" = names(pk_out),
                 "param_value" = pk_out)

  outval


}

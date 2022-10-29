#' Assign an estimated limit of quantification wherever missing
#'
#' Sometimes in mass spectrometry a compound is detected, but the signal is
#' too low and falls outside the range that can be interpreted quantitatively
#' (that is, as a concentration). These detects can still be used to refine
#' parameter estimates by preferring curve fits that predict concentrations
#' below the known "limit of quantitation" (LOQ) -- the lowest detectable quantifiable
#' concentration. LOQ's are not always available from data reported in the
#' literature, so this function estimates a LOQ such that all the reported
#' values are above the LOQ. The LOQ is assumed to vary from study to study,
#' as each study would have different mass spectrometry equipment and methods.
#' LOQ also is assumed to vary by media, as preparation of samples may be different.
#'
#' @param cvt.data A table of concentration vs. time data with at least the
#' following columns:
#'
#' @param source.col The column indicating different research studies (and
#' there for different mass spectrometry methods)
#'
#' @param chem.col The column indicating chemical identity
#'
#' @param media.col The column indicating sample medium
#'
#' @param value.col The column indicating concentration
#'
#' @param calc.loq.col The column to be added to the table holding the estimated
#' LOQ's
#'
#' @return A table with column "calc_loq" added
#'
#' @author John Wambaugh
#'
#' @export estimate_loq
estimate_loq <- function(cvt.data,
                         source.col = "fk_study_id",
                         chem.col = "dsstox_casrn",
                         media.col = "conc_medium_normalized",
                         value.col = "conc",
                         calc.loq.col = "calc_loq")
{

  for (this.study in unique(cvt.data[,source.col]))
  {
    this.subset1 <- subset(cvt.data,cvt.data[,source.col]==this.study)
    for (this.chem in unique(this.subset1[,chem.col]))
    {
      this.subset2 <- subset(this.subset1,this.subset1[,chem.col]==this.chem)
      for (this.media in unique(this.subset2[,media.col]))
      {
        this.subset3 <- subset(this.subset2,this.subset2[media.col]==this.media &
                                 !is.na(this.subset2[,value.col]))
        this.subset3 <- subset(this.subset3,this.subset3[,value.col]>0)
        this.loq <- suppressWarnings(min(this.subset3[,value.col], na.rm=TRUE))
        this.loq <- this.loq*0.9
        cvt.data[
          cvt.data[,source.col]==this.study &
            cvt.data[,chem.col]==this.chem &
            cvt.data[,media.col]==this.media,calc.loq.col] <- this.loq
      }
    }
  }
  cvt.data <- subset(cvt.data,!is.infinite(cvt.data[,calc.loq.col]))

  return(cvt.data)
}

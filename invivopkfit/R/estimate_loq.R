#' Assign an estimated limit of quantification wherever missing
#'
#' Sometimes in mass spectrometry a compound is detected, but the signal is too
#' low and falls outside the range that can be interpreted quantitatively (that
#' is, as a concentration). These detects can still be used to refine parameter
#' estimates by preferring curve fits that predict concentrations below the
#' known "limit of quantitation" (LOQ) -- the lowest detectable quantifiable
#' concentration. LOQ's are not always available from data reported in the
#' literature, so this function estimates a LOQ such that all the reported
#' values are above the LOQ.

#' The LOQ is assumed to vary from study to study, as each study would have
#' different mass spectrometry equipment and methods. LOQ also is assumed to
#' vary by media, as preparation of samples may be different. For each unique
#' combination of study, chemical, and medium, the LOQ is estimated by
#' multiplying the lowest detected concentration (the minimum non-NA reported
#' concentration greater than zero) by a constant factor (\code{loq.factor}). If
#' there are no detected concentrations for a given combination of study,
#' chemical, and medium (i.e., all reported concentrations are NA or 0), then
#' that study/chemical/medium combination will be removed from the returned
#' data.frame.
#'
#' @param cvt.data A \code{data.frame} of concentration vs. time data with at
#'   least the following columns:
#'
#' @param source.col The name of the column indicating different research study
#'   source documents (and therefore different mass spectrometry methods).
#'   Default "document_id".
#'
#' @param chem.col The name of the column indicating chemical identity. Default
#'   "dsstox_casrn".
#'
#' @param media.col The name of the column indicating sample medium. Default
#'   "conc_medium_normalized".
#'
#' @param value.col The name of the column indicating concentration. Default
#'   "conc". Nondetects should be reported as 0 or NA in this column.
#'
#' @param calc.loq.col The name of the column to be added to the table holding
#'   the estimated LOQ's. Default "calc_loq".
#'
#' @param loq.factor A factor by which to multiply the lowest detected value to
#'   estimate an LOQ. Default 0.45.
#'
#' @return A \code{data.frame} that is the same as \code{cvt.data}, but with the
#'   column specifed in \code{calc_loq} added, and any rows for
#'   source/chemical/medium combinations with no detected concentrations removed
#'
#' @author John Wambaugh
#'
#' @export estimate_loq
estimate_loq <- function(cvt.data,
                         source.col = "document_id",
                         chem.col = "dsstox_casrn",
                         media.col = "conc_medium_normalized",
                         value.col = "conc",
                         calc.loq.col = "calc_loq",
                         loq.factor = 0.45)
{

  for (this.study in unique(cvt.data[,source.col]))
  {
    this.subset1 <- subset(cvt.data,
                           cvt.data[,source.col]==this.study)
    for (this.chem in unique(this.subset1[,chem.col]))
    {
      this.subset2 <- subset(this.subset1,
                             this.subset1[,chem.col]==this.chem)
      for (this.media in unique(this.subset2[,media.col]))
      {
        this.subset3 <- subset(this.subset2,
                               this.subset2[media.col]==this.media &
                                 !is.na(this.subset2[,value.col])
                               )
        this.subset3 <- subset(this.subset3,
                               this.subset3[,value.col]>0)
        this.loq <- suppressWarnings(min(this.subset3[,value.col],
                                         na.rm=TRUE))
        this.loq <- this.loq*loq.factor
        cvt.data[
          cvt.data[,source.col]==this.study &
            cvt.data[,chem.col]==this.chem &
            cvt.data[,media.col]==this.media,calc.loq.col] <- this.loq
      }
    }
  }
  cvt.data <- subset(cvt.data,
                     !is.infinite(cvt.data[,calc.loq.col])
                     )

  return(cvt.data)
}

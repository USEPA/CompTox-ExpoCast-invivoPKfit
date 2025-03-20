# invivoPKfit 2.0.1

* Ensured only optimized parameters will be used when calculating AIC. This is really
only a concern when using models that hold many of the total parameters constant.  
* Added a new feature allowing users to calculate log-likelihood from predictions
generated outside invivoPKfit.
* Updated documentation and cleaned up vignettes.


# invivoPKfit 2.0.0

* Initial CRAN submission.  
* Redesigned package to an S3 object-oriented approach with some ggplot2-like syntax.  
* Re-factored function/method code to ensure the code base is more legible to future contributors.  
* Implemented messages that save and update data processing "status".  
* Parallel processing parameter optimization was added as an option to take advantage of multi-processing cores.  
* Included more chemical records than previously published in the `cvt` object, represented a processed portion of the CvTdb database.

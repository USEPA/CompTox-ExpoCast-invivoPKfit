# invivoPKfit 2.0.2

* Minor bug fix for adjusting rate constants after fitting to standard "1/hour".  
* Improved `pk` object organization and CLI message outputs.  
* Externalized parallel processing for fitting with `mirai` and `purrr` 1.1.0.
Now, users should set the number of processes through those interfaces if
parallel or distributed computing is desired.  
* Setting data, error, NCA, sd, and LOQ groups now uses unquoted data variable names,
the latter two have their own separate `pkproto` object setups: `stat_sd_group()` and `stat_loq_group()`.  
* Parameter starts and and estimates are stored as singular values in a `list-column`,
which allows character vector constant parameters for calls to more complex PBTK models,
see `model_httk_gas_pbtk`.  
* Users may now set specific parameters to optimize and specific start values using the new functions
`set_params_optimize()` and `set_params_starts()`.
* Updated `cvt` object to align with latest CvTdb release.


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

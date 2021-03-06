% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analytic_model_fun.R
\name{analytic_model_fun}
\alias{analytic_model_fun}
\title{Evaluate analytic models}
\usage{
analytic_model_fun(params, dose, times, time.units, iv.dose, model)
}
\arguments{
\item{params}{A named list of parameter values. Must match the parameters of \code{model}.}

\item{dose}{A dose in units of mg/kg.}

\item{times}{A vector of times in hours or days. Ideally should be sorted.}

\item{time.units}{The units of \code{times}: "h" for hours, "d" for days.}

\item{iv.dose}{TRUE for IV dosing, FALSE for PO dosing.}

\item{model}{Analytic model to evaluate. Currently only "1compartment" or
"2compartment" are implemented.}
}
\value{
A matrix with three columns:
\describe{
\item{time}{Time in days}
\item{Ccompartment}{Plasma concentrations}
\item{AUC}{Area under the curve}}
}
\description{
Evaluate analytic models
}
\author{
Caroline Ring
}

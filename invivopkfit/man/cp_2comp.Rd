% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cp_2comp.R
\name{cp_2comp}
\alias{cp_2comp}
\title{Analytical 2-compartment model}
\usage{
cp_2comp(params, time, dose, iv.dose)
}
\arguments{
\item{params}{A named list of parameter values including the following: 
  \describe{ \item{k12}{Rate at which compound moves from central to 
  peripheral compartment} \item{k21}{Rate at which compound moves from 
  peripheral to central compartment} \item{kelim}{Elimination rate} 
  \item{V1}{Apparent volume of central compartment} } For oral administration
  (\code{iv.dose} FALSE), \code{params} must also include: \describe{ 
  \item{Fgutabs}{Oral bioavailability} \item{kgutabs}{Rate of absorption from
  gut} }

author Caroline Ring, John Wambaugh}

\item{time}{A vector of time values, in hours}

\item{dose}{A dose in mg/kg}

\item{iv.dose}{TRUE for single IV bolus dose, FALSE for single oral dose}
}
\value{
A vector of plasma concentration values corresponding to each value 
  in \code{time}
}
\description{
Analytical 2-compartment model
}

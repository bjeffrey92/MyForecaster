% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/forecast_comparisons.R
\name{check_stationary}
\alias{check_stationary}
\title{Checks if data is stationary}
\usage{
check_stationary(time_series)
}
\arguments{
\item{time_series}{any time series object}
}
\description{
uses an augmented dicky fuller test to determine whether the time series is stationary, if not differences until it is
}

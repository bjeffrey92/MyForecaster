% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/forecast_comparisons.R
\name{PI_accuracy}
\alias{PI_accuracy}
\title{Accuracy of model produced prediction intervaks}
\usage{
PI_accuracy(testing_data, lower, upper)
}
\arguments{
\item{testing_data}{the end of the time series that was removed for forecasting}

\item{lower}{a lower prediction interval}

\item{upper}{an upper prediction interval}
}
\description{
returns the percentage of data points in testing data between lower and upper PI
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/forecasting_functions.R
\name{ets_forecast}
\alias{ets_forecast}
\title{Exponential Smoothing Forecast}
\usage{
ets_forecast(training_data, forecasting_horizon, frequency, PI = 95,
  apply_limits = TRUE)
}
\arguments{
\item{training_data}{time series to be forecasted}

\item{forecasting_horizon}{time in units of training_data to be forecasted}

\item{frequency}{12 for months, 4 for quarters etc}

\item{PI}{width of prediction intervals, default is 95 (ie 95 percent)}

\item{apply_limits}{keep data with 0 and 100, default is TRUE}
}
\description{
uses AICc to find best fitting ets model and then applies this to forecast function
}

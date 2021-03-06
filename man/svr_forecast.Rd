% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/forecasting_functions.R
\name{svr_forecast}
\alias{svr_forecast}
\title{Support Vector Regression Forecast}
\usage{
svr_forecast(training_data, forecasting_horizon, frequency,
  apply_limits = TRUE)
}
\arguments{
\item{training_data}{time series to be forecasted}

\item{forecasting_horizon}{time in units of training_data to be forecasted}

\item{frequency}{12 for months, 4 for quarters etc}

\item{apply_limits}{keep data with 0 and 100, default is TRUE}
}
\description{
fits an SVR model and extrapolates from it to produce forecasts
}

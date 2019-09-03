library(forecast)
library(ggplot2)
library(stringr)
library(TSPred)
library(tseries)

#' My adaptation of check_residuals function
#'
#' @description adapted from forecast::checkresiduals
#' @export
#' @param object a fitted forecasting model

my_checkresiduals <- function (object, lag, df = NULL, test, plot = TRUE, ...){

    showtest <- TRUE
    if (missing(test)) {
        if (is.element("lm", class(object))) {
            test <- "BG"
        }
        else {
            test <- "LB"
        }
        showtest <- TRUE
    }
    else if (test != FALSE) {
        test <- match.arg(test, c("LB", "BG"))
        showtest <- TRUE
    }
    else {
        showtest <- FALSE
    }
    if (is.element("ts", class(object)) | is.element("numeric",
        class(object))) {
        residuals <- object
        object <- list(method = "Missing")
    }
    else {
        residuals <- residuals(object)
    }
    if (length(residuals) == 0L) {
        stop("No residuals found")
    }
    if ("ar" %in% class(object)) {
        method <- paste("AR(", object$order, ")", sep = "")
    }
    else if (!is.null(object$method)) {
        method <- object$method
    }
    else if ("HoltWinters" %in% class(object)) {
        method <- "HoltWinters"
    }
    else if ("StructTS" %in% class(object)) {
        method <- "StructTS"
    }
    else {
        method <- try(as.character(object), silent = TRUE)
        if ("try-error" %in% class(method)) {
            method <- "Missing"
        }
        else if (length(method) > 1 | base::nchar(method[1]) >
            50) {
            method <- "Missing"
        }
    }
    if (method == "Missing") {
        main <- "Residuals"
    }
    else {
        main <- paste("Residuals from", method)
    }
    if (plot) {
        suppressWarnings(ggtsdisplay(residuals, plot.type = "histogram",
            main = main, ...))
    }
    if (is.element("forecast", class(object))) {
        object <- object$model
    }
    if (is.null(object) | !showtest) {
        return(invisible())
    }
    freq <- frequency(residuals)
    if (is.element("ets", class(object))) {
        df <- length(object$par)
    }
    else if (is.element("Arima", class(object))) {
        df <- length(object$coef)
    }
    else if (is.element("bats", class(object))) {
        df <- length(object$parameters$vect) + NROW(object$seed.states)
    }
    else if (is.element("lm", class(object))) {
        df <- length(object$coefficients)
    }
    else if (method == "Mean") {
        df <- 1
    }
    else if (grepl("Naive", method, ignore.case = TRUE)) {
        df <- 0
    }
    else if (method == "Random walk") {
        df <- 0
    }
    else if (method == "Random walk with drift") {
        df <- 1
    }
    else {
        df <- NULL
    }
    if (missing(lag)) {
        lag <- ifelse(freq > 1, 2 * freq, 10)
        lag <- min(lag, length(residuals)/5)
        lag <- max(df + 3, lag)
    }
    if (!is.null(df)) {
        if (test == "BG") {
            BGtest <- lmtest::bgtest(object, order = lag)
            BGtest$data.name <- main
            return(BGtest)
        }
        else {
            LBtest <- Box.test(zoo::na.approx(residuals), fitdf = df,
                lag = lag, type = "Ljung")
            LBtest$method <- "Ljung-Box test"
            LBtest$data.name <- main
            names(LBtest$statistic) <- "Q*"
            return(LBtest)
        }
    }
}

#' My adaption of autolayer function
#'
#' @description adapted from forecast:::autolayer to allow objects of length 1 to be plotted also
#' @export
#' @param object a mean of a forecast object
#' @param series name of series to be plotted

my_autolayer <- function(object, series = NULL, ...){
  
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("ggplot2 is needed for this function to work. Install it via install.packages(\"ggplot2\")",
          call. = FALSE)
  }
  else {
      tsdata <- data.frame(timeVal = as.numeric(time(object)),
          series = ifelse(is.null(series), deparse(substitute(object)),
              series), seriesVal = as.numeric(object))
      if (length(object) > 1) {
          ggplot2::geom_line(ggplot2::aes_(x = ~timeVal, y = ~seriesVal,
              group = ~series, colour = ~series), data = tsdata,
              ...)
      }
      else {
          ggplot2::geom_point(ggplot2::aes_(x = ~timeVal, y = ~seriesVal,
            group = ~series, colour = ~series), data = tsdata, ...)
      }
  }
}


#' Mean error of forecast
#'
#' @description calculates mean error of a forecast
#' @export
#' @param testing_data the end of the time series that was removed for forecasting
#' @param point_forecasts mean from a forecast object

ME <- function(testing_data, point_forecasts){

   mean_error <- mean(testing_data - as.numeric(point_forecasts))
   mean_error <- round(mean_error, digits = 3)
   return(mean_error)
}


#' Symmetirc mean absolute percentage error
#'
#' @description calculates the symmetric mean absolute percentage error
#' @export
#' @param testing_data the end of the time series that was removed for forecasting
#' @param point_forecasts mean from a forecast object

SMAPE <- function(testing_data, point_forecasts){
  
  SMAPE <- TSPred::sMAPE(testing_data, point_forecasts)
  return(SMAPE)
}


#' Root Mean Squared Error
#'
#' @description calculates root mean squared error of a forecast
#' @export
#' @param testing_data the end of the time series that was removed for forecasting
#' @param point_forecasts mean from a forecast object

RMSE <- function(testing_data, point_forecasts){
  
  root_mean_squared_error <- sqrt(mean((testing_data -
                                  as.numeric(point_forecasts)) ^ 2))
  root_mean_squared_error <- round(root_mean_squared_error, digits = 3)
  return(root_mean_squared_error)
}


#' Prediction Interval Width
#'
#' @description calculates the width of the forecast prediction intervals
#' @export
#' @param lower a lower prediction interval
#' @param upper an upper prediction interval

PI_width <- function(lower, upper){

  mean_width <- mean(upper - lower)
  return(mean_width)
}


#' Accuracy of model produced prediction intervaks
#'
#' @description returns the percentage of data points in testing data between lower and upper PI
#' @export
#' @param testing_data the end of the time series that was removed for forecasting
#' @param lower a lower prediction interval
#' @param upper an upper prediction interval

PI_accuracy <- function(testing_data, lower, upper){

  x <- testing_data[testing_data >= lower & testing_data <= upper]
  percent_accurate <- length(x)/length(testing_data) * 100
  return(percent_accurate)
}


#' Inspects Residuals of Model fit
#'
#' @description returns the mean of the residuals and a p value of pormanteau test for normallity
#' @export
#' @param fit a fitted forecsting model
#' @param residuals the residuals of the fitted model
#' @param forecast_name name of the forecast object

inspect_residuals <- function(fit, residuals, forecast_name){

  residuals_mean <- mean(residuals, na.rm = TRUE)

  #for these forecasting algorithms checkresiduals func will interpret fit to estimate df
  methods_1 <- c('mean_forecast',
                'naive_forecast',
                'naive_forecast_with_drift',
                'ets_forecast',
                'arima_forecast')
  if (forecast_name %in% methods_1){
    p_value <- my_checkresiduals(fit, plot = FALSE)$p.value
  } else{
    p_value <- NA
  }

  output_list <- c(residuals_mean, p_value)
  return(output_list)
}

#' Checks if data is stationary
#'
#' @description uses an augmented dicky fuller test to determine whether the time series is stationary, if not differences until it is
#' @export
#' @param time_series any time series object 

check_stationary <- function(time_series){

  p <- tseries::adf.test(time_series)$p.value
  differencing_order <- 0
  while (p > 0.05){
    time_series <- diff(time_series)
    p <- tseries::adf.test(time_series)$p.value
    differencing_order <- differencing_order + 1
  }
  output <- list(time_series, differencing_order)
  return(output)
}

#' Reverse differencing
#'
#' @description reverts differencing that was applied to make training data stationary
#' @export
#' @param time_series any time series object 
#' @param differencing_order order of differencing that was applied

reverse_differencing <- function(time_series, differencing_order){

  output <- diffinv(time_series, differences=differencing_order)
  return(output)
}


#' Compare a number of forecasting methods on one time series
#'
#' @description passes time series produced by convert_to_raw_time_series to different forecasting algorithms and then collates the results
#' @export
#' @param all_forecasting_functions list of lists describing the names of forecasting functions to be compared. Follow structure of list of same name exported from module 
#' @param time_series output of convert_to_raw_time_series
#' @param location where was data taken from for plot title
#' @param ab name of antibiotic for plot title
#' @param organism name of organism for plot title
#' @param frequency 12 for months, 4 for quarters
#' @param apply_limits true or false, limit forecasts to 0-100

do_forecasts <- function(all_forecasting_functions, time_series, forecasting_horizon,
                        max_forecasting_horizon, location, ab, organism,
                        frequency, apply_limits = TRUE){

  last_date <- time(time_series)[[length(time_series)]]
  training_data <- stats::window(time_series, end = last_date -
                                max_forecasting_horizon) #splits data for testing
  testing_data <- stats::window(time_series, start = last_date -
                                max_forecasting_horizon + 1/frequency,
                                end = last_date - max_forecasting_horizon +
                                forecasting_horizon)

  LB_p_value <- Box.test(training_data, type = 'Lj')$p.value #Ljung-Box test to assess difference from white noise
  LB_p_value <- round(LB_p_value, digits = 3)

  plot <- forecast::autoplot(time_series) #base layer of plot to be added to

  output_list <- list(organism,
                      location,
                      ab,
                      forecasting_horizon,
                      LB_p_value,
                      round(mean(training_data), digits = 3),
                      round(mean(diff(training_data)), digits = 3),
                      round(var(training_data), digits = 3))

  ##APPLY ALTERNATIVE FORECASTING ALGORITHMS AND CALCULATE THEIR RMSE AND ME
  for (forecast_name in names(all_forecasting_functions)){

    #foreacst will be list where first item is point foreast followed by lower and then upper 95% prediction intervals
    forecast_func <- all_forecasting_functions[[forecast_name]][[1]]
    must_be_stationary <- all_forecasting_functions[[forecast_name]][[2]]

    if (must_be_stationary){
      td <- check_stationary(training_data)
      training_data <- td[[1]]
      differencing_order <- td[[2]]
    } else{
      differencing_order <- NA
    }

    forecast <- get(forecast_func)(training_data,
      forecasting_horizon * frequency, frequency, apply_limits) #get evaluates function text stored in functions list

    if (!is.na(differencing_order)){
      for (i in 1:length(forecast)){
        forecast[[i]] <- reverse_differencing(foreacst[[i]], differencing_order)
      }
    }

    print(forecast_name)
    if (!is.na(forecast[[1]][[1]])){
      forecast_ME <- ME(testing_data, forecast[[1]]) #calculate mean error
      forecast_RMSE <- RMSE(testing_data, forecast[[1]]) #calculate root mean squared error
      forecast_PI_Width <- PI_width(forecast[[2]], forecast[[3]])
      forecast_PI_Accuracy <- PI_accuracy(testing_data, forecast[[2]],
                                          forecast[[3]])
      forecast_sMAPE <- SMAPE(testing_data, forecast[[1]])
      residuals_analysis <- inspect_residuals(forecast[[4]], forecast[[5]],
                                              forecast_name)
      residuals_mean <- residuals_analysis[[1]]
      pormanteau_test_p_value <- residuals_analysis[[2]]

      plot <- plot + my_autolayer(forecast[[1]], series = forecast_name) #add layer to plot
    } else{
      forecast_ME <- NA
      forecast_RMSE <- NA
      forecast_PI_Width <- NA
      forecast_PI_Accuracy <- NA
      forecast_sMAPE <- NA
      residuals_mean <- NA
      pormanteau_test_p_value <- NA
    }
    output_list[[paste(forecast_name, 'ME', sep = '_')]] <- forecast_ME
    output_list[[paste(forecast_name, 'RMSE', sep = '_')]] <- forecast_RMSE
    output_list[[paste(forecast_name, 'sMAPE', sep = '_')]] <- forecast_sMAPE
    output_list[[paste(forecast_name, 'PI_width', sep = '_')]] <- forecast_PI_Width
    output_list[[paste(forecast_name, 'PI_accuracy', sep = '_')]] <- forecast_PI_Accuracy
    output_list[[paste(forecast_name, 'residuals_mean', sep = '_')]] <- residuals_mean
    output_list[[paste(forecast_name, 'pormanteau_test_p_value', sep = '_')]] <- pormanteau_test_p_value
  }

  ##PLOT FORECASTS AND REAL DATA
  plot_title <- paste(organism, 'percentage resistance to', ab, 'in', location,
                      '(h =', forecasting_horizon, ')')
  filename <- paste0(organism, '_', location, '_', ab, '_forecasting_horizon',
                      forecasting_horizon, '.png')
  filename <- gsub(' ', '_', filename) #remove spaces from file name
  plot <- plot + ggplot2::xlab('Year') + ggplot2::ylab('Percetage Resistance') +
          ggplot2::ggtitle(plot_title)
  ggplot2::ggsave(filename)

  return(output_list)
}
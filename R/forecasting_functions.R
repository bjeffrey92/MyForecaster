library(forecast)
library(forecTheta)
library(EMD)
library(fractaldim)
library(plyr)
library(dplyr)
library(e1071)
library(prophet)

#' Limit output to between 0 and 100 for data which is a percentage 
#'
#' @description ensures that predicted values are between 0 and 100
#' @export 
#' @param forecast a forecast object

apply_limits <- function(forecast){

  for (number in 1:5){
    if (number != 4){
      item <- forecast[[number]]
      for (horizon in 1:length(item)){
        horizon_value <- item[[horizon]]
        if (is.na(horizon_value)){  
          #do nothing 
        } else if (horizon_value < 0){
          horizon_value <- 0
        } else if (horizon_value > 100){
          horizon_value <- 100
        }
        item[[horizon]] <- horizon_value
      }
      forecast[[number]] <- item
    }
  }
  return(forecast)
}

#' Mean Forecast
#'
#' @description mean of all previous points
#' @export
#' @param training_data time series to be forecasted
#' @param forecasting_horizon time in units of training_data to be forecasted
#' @param frequency 12 for months, 4 for quarters etc
#' @param apply_limits keep data with 0 and 100, default is TRUE 

mean_forecast <- function(training_data, forecasting_horizon, frequency,
                          apply_limits = TRUE){

  mean_forecast <- forecast::meanf(training_data, forecasting_horizon)
  output <- list(mean_forecast[]$mean, mean_forecast[]$lower[,2], #returns point forecast and upper and lower 95% ;rediction intervals
                mean_forecast[]$upper[,2], mean_forecast,
                mean_forecast[]$residuals)

  if (apply_limits){
    output <- apply_limits(output)
  }
  return(output)
}

#' Naive Forecast
#'
#' @description all equal to last point
#' @export
#' @inheritParams mean_forecast


naive_forecast <- function(training_data, forecasting_horizon, frequency,
                            apply_limits = TRUE){

  naive_forecast <- forecast::rwf(training_data, forecasting_horizon, drift = FALSE)
  output <- list(naive_forecast[]$mean, naive_forecast[]$lower[,2], #returns point forecast and upper and lower 95% ;rediction intervals
                naive_forecast[]$upper[,2], naive_forecast,
                naive_forecast[]$residuals)

  if (apply_limits){
    output <- apply_limits(output)
  }
  return(output)
}
  
#' Naive Forecast With drift
#'
#' @description same as naive but drift component is fitted
#' @export
#' @inheritParams mean_forecast

naive_forecast_with_drift <- function(training_data, forecasting_horizon,
                                      frequency, apply_limits = TRUE){

  if (var(training_data) == 0){ #will not work if all values in training set are equal
    return(NA)
  } else {
    naive_forecast_with_drift <- forecast::rwf(training_data, forecasting_horizon,
      drift = TRUE)
    output <- list(naive_forecast_with_drift[]$mean,
                  naive_forecast_with_drift[]$lower[,2], #returns point forecast and upper and lower 95% ;rediction intervals
                  naive_forecast_with_drift[]$upper[,2],
                  naive_forecast_with_drift,
                  naive_forecast_with_drift[]$residuals)

    if (apply_limits){
      output <- apply_limits(output)
    }
    return(output)
  }
}



#' Exponential Smoothing Forecast
#'
#' @description uses AICc to find best fitting ets model and then applies this to forecast function
#' @export
#' @inheritParams mean_forecast

ets_forecast <- function(training_data, forecasting_horizon, frequency,
                        apply_limits = TRUE){

  fit <- forecast::ets(training_data)
  ets_forecast <- forecast::forecast(fit, h = forecasting_horizon)
  output <- list(ets_forecast[]$mean, ets_forecast[]$lower[,2], #returns point forecast and upper and lower 95% ;rediction intervals
                ets_forecast[]$upper[,2], fit, fit$residuals)

  if (apply_limits){
    output <- apply_limits(output)
  }
  return(output)
}



#' ARIMA Forecast
#'
#' @description uses AICc to find best fitting arima model and then applies this to forecast function
#' @export
#' @inheritParams mean_forecast

arima_forecast <- function(training_data, forecasting_horizon, frequency,
                          apply_limits = TRUE){

  fit <- forecast::auto.arima(training_data)
  arima_forecast <- forecast::forecast(fit, h = forecasting_horizon)
  output <- list(arima_forecast[]$mean, arima_forecast[]$lower[,2], #returns point forecast and upper and lower 95% ;rediction intervals
                arima_forecast[]$upper[,2], fit, fit$residuals)

  if (apply_limits){
    output <- apply_limits(output)
  }
  return(output)
}



#' Neural Network Autoregression Forecast
#'
#' @description fits an autoregressive feed forward neural network model to the data and estimates uncertainity with bootstrapping
#' @export
#' @inheritParams mean_forecast

nn_autoregression_forecast <- function(training_data, forecasting_horizon,
                                      frequency, apply_limits = TRUE){


  fit <- forecast::nnetar(training_data)
  nnetar_forecast <- forecast::forecast(fit, PI = TRUE, h = forecasting_horizon) #fits model and generates a point forecast

  output <- list(nnetar_forecast[]$mean, nnetar_forecast[]$lower[,2], #returns point forecast and upper and lower 95% ;rediction intervals
                nnetar_forecast[]$upper[,2], fit, fit$residuals)

  if (apply_limits){
    output <- apply_limits(output)
  }
  return(output)
  #implemented within the forecast function
  # sim <- ts(matrix(0, nrow=forecasting_horizon, ncol=30), #empty matrix to populate with bootstrap samples to generate an approximate prediction interval
  #         start=end(training_data)[1]+1)
  # for (i in seq(30)){
  #   sim[,i] <- simulate(fit, nsim = forecasting_horizon) #does 30 forward projections to estimate uncertainty
  # }
}


#' Support Vector Regression Forecast
#'
#' @description fits an SVR model and extrapolates from it to produce forecasts
#' @export
#' @inheritParams mean_forecast

svr_forecast <- function(training_data, forecasting_horizon, frequency,
                        apply_limits = TRUE){

  #formats the data
  data <- cbind(as.vector(training_data), 1:length(training_data))
  data <- as.data.frame(data)
  names(data) <- c('Y', 'X')

  #performs grid search to optimise hyper parameters
  model <- e1071::tune(svm, Y ~ X,  data = data,
              ranges = list(epsilon = seq(0,1,0.1), cost = 2^(2:9)))

  #need df to fill in for forward predictions
  last_training_point <- tail(data, n = 1)$X
  prediction_df <- cbind(1:length(forecasting_horizon),
                        (last_training_point + 1):(last_training_point +
                          forecasting_horizon))
  prediction_df <- as.data.frame(prediction_df)
  names(prediction_df) <- c('Y', 'X')

  #extrapolates from regression function to predict future values
  prediction <- predict(model$best.model, prediction_df)
  return(prediction)
}



#' Standard Theta Forecast
#'
#' @description fits a theta model first described by V.Assimakopoulos annd K.Nikolopoulos, 2000
#' @export
#' @inheritParams mean_forecast

standard_theta_forecast <- function(training_data, forecasting_horizon,
                                    frequency, apply_limits = TRUE){

  if (forecasting_horizon == 1){ #bug in package means can't do forecast with h = 1
  last_date <- time(training_data)[[length(training_data)]]
  forecast_date <- last_date + 1

  forecasting_horizon <-2
  theta_forecast <- forecTheta::stm(training_data, h = forecasting_horizon)

  mean <- ts(theta_forecast$mean[1], start = forecast_date)
  lower <- ts(theta_forecast$lower[,3][1], start = forecast_date)
  upper <- ts(theta_forecast$upper[,3][1], start = forecast_date)

  output <- list(mean, lower, upper)
  } else{
    theta_forecast <- forecTheta::stm(training_data, h = forecasting_horizon)
    output <- list(theta_forecast$mean, theta_forecast$lower[,3],
                  theta_forecast$upper[,3], theta_forecast,
                  theta_forecast$residuals)
  }

  if (apply_limits){
    output <- apply_limits(output)
  }
  return(output)
}


#' Dynamic Optimised Theta Forecast
#'
#' @description fits dynamic optimised theta model as developed by Fiorucci et al 2016
#' @export
#' @inheritParams mean_forecast

dotm_forecast <- function(training_data, forecasting_horizon, frequency,
                          apply_limits = TRUE){

  if (forecasting_horizon == 1){ #bug in package means can't do forecast with h = 1
    last_date <- time(training_data)[[length(training_data)]]
    forecast_date <- last_date + 1

    forecasting_horizon <-2
    theta_forecast <- forecTheta::dotm(training_data, h = forecasting_horizon)

    mean <- ts(theta_forecast$mean[1], start = forecast_date)
    lower <- ts(theta_forecast$lower[,3][1], start = forecast_date)
    upper <- ts(theta_forecast$upper[,3][1], start = forecast_date)

    output <- list(mean, lower, upper)
  } else{
    theta_forecast <- forecTheta::dotm(training_data, h = forecasting_horizon)
    output <- list(theta_forecast$mean, theta_forecast$lower[,3],
                  theta_forecast$upper[,3], theta_forecast,
                  theta_forecast$residuals)
  }

  if (apply_limits){
    output <- apply_limits(output)
  }
  return(output)
}


#' Teodoro and Lovis Forecast
#'
#' @description fits the model developed by Teodoro and Lovis, 2013
#' @export
#' @inheritParams mean_forecast
#' @param filter filter to be used for selecting components of EMD for KNN forecasting. Default is DECF

t_and_l_forecast <- function(training_data, forecasting_horizon, frequency,
                            apply_limits = TRUE, filter = 'DECF'){

  filter_types <- c('DECA','DECF', 'DECS')
  if (!(filter %in% filter_types)){
    filter <- as.character(filter)
    message <- paste('ERROR, unknown filter type:', filter)
    stop(message)
  }

  ################
  #MODEL TRAINING#
  ################
  #decomposes the time series and estimates m and k

  #empirical mode decomposition of time series
  decomposed_ts <- EMD::emd(training_data)# performs emd, returns all IMFs and the residual
  components <- list()

  #if there are no imfs then the residue will be taken as the only component for forecasting
  if (decomposed_ts$nimf == 0){
    components[[1]] <- decomposed_ts$residue
  } else{
    for (i in 1:ncol(decomposed_ts$imf)){
      components[[length(components) + 1]] <- decomposed_ts$imf[,i] #appends each imf to components list
    }
    components[[length(components) + 1]] <- decomposed_ts$residue
  }

  #model determines which components are included, this is called filter in t and l paper
  count_zero_crossings <- function(signal){
    zero_crossings <- 0
    for (i in 2:length(signal)){
      if (abs(signal[i] - signal[i-1]) > abs(signal[i])){
        zero_crossings <- zero_crossings + 1
      }
    }
    return(zero_crossings)
  }

  if (filter == 'DECF'){ #DECF filters out low frequency components defined as those with period less than 10 weeks or 2.5 months
    filtered_components <- list()
    for (component in components){
      include <- FALSE
      zero_crossings <- count_zero_crossings(component)
      if (frequency == 12){ #if data is monthly
        if (zero_crossings < length(component)/2.5){
          include <- TRUE
        }
      } else if (frequency == 52){ #if data is weekly
        if (zero_crossings < length(component)/10){
          include <- TRUE
        }
      } else if (frequency == 4){ #if data is quarterly
        if (zero_crossings < length(component)){ #approximate to ten week filter
          include <- TRUE
        }
      } else{
        stop('Unknown frequncy value for DECF filter')
      }
      if (include){
        filtered_components[[length(filtered_components) + 1]] <- component
      }
    }
    components <- filtered_components
  }
  # if (filter == 'DECS'){ #use statistical test derived by Wu and Huang to decide which components to include
  #   stop('DECS filter is not implemented')
  #   if (length(components) > 2){
  #     stop('Cannot use DECS filter, only one IMF detected')
  #   }
  #   filtered_components <- list()
  #   first_component <- components[[1]] #assume first component is noise
  #
  #   T_n <- 1/count_zero_crossings(first_component)
  #   sigma <- (2*T_n)/length(first_component)
  #   upper_spread <- -log(T_n) + 2*sigma
  #   lower_spread <- -log(T_n) - 2*sigma
  #
  #   for (component in tail(components, length(components) - 1)){
  #     include <- FALSE
  #     T_n <- 1/count_zero_crossings(component)
  #     sigma <- (2*T_n)/length(component)
  #     log_10_variance <- log10(sigma^2)
  #     log_10_period <- log10(T_n)
  #   }
  # }

  #following method in http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.68.1974&rep=rep1&type=pdf
  #estimates fractal dim at each lag 1-10, where slope levels (defined in above paper) take as optimal lag length

  m_max <- 10 #max dimension
  component_optimal_m_values <- c()
  for (component in components){
    all_lagged_components <- list()
    all_lagged_components[[1]] <- component

    lagged_componnent_fl_values <- c()
    for (m in 1:m_max){

      lagged_data <- dplyr::lag(component, m) #get m step lag
      all_lagged_components[[m + 1]] <- lagged_data
      lag_matrix <- do.call('cbind', all_lagged_components) #create matrix of all lagged data
      lag_matrix <- na.omit(lag_matrix) #drop rows with na

      f_l <- fractaldim::fd.estim.boxcount(lag_matrix)$fd #calculate fractal dimension of matrix using box counting algorithm
      lagged_componnent_fl_values <- c(lagged_componnent_fl_values, f_l)
    }

    #estimate point where slope flattens using definition in paper cited above
    #point where difference in consecutive values equals max(0.3, 10% of ma of fractal dim)
    found_optimal_m <- FALSE
    for (i in 1:length(lagged_componnent_fl_values)){
      if (i > 1){
        dif <- abs(lagged_componnent_fl_values[i] -
                  lagged_componnent_fl_values[i-1]) #diff between two consecutive values
        MA <- mean(c(lagged_componnent_fl_values[i],
                      lagged_componnent_fl_values[i-1])) #moving average of fd

        if (dif <= MA || dif <= 0.3){ #has slope levelled off
          component_optimal_m <- i
          component_optimal_m_values <- c(component_optimal_m_values,
                                          component_optimal_m)
          found_optimal_m <- TRUE
          break
        }
      }
    }
    if (!found_optimal_m){ #if slope doesn't level off take optimal m as m_max
      component_optimal_m_values <- c(component_optimal_m_values, m_max)
    }
  }

  m <- mean(component_optimal_m_values) #optimal value of m

  #compute optimal k using cross validation
  for (j in m:(length(training_data) -1)){
    all_optimal_k_values <- c()
    for (component in components){
      if (j > m){
        #delay vector and next point in component of time series
        delay_vector <- c()
        for (i in 1:m){
          delay_vector <- c(delay_vector, component[[j - i]])
        }
        x_ij_plus_1 <- component[[j]]
        all_distances <- list()

        for (k in 1:length(delay_vector)){ #compute euclidean distance from target point for each set of knn
          distance <- abs(x_ij_plus_1 - mean(head(delay_vector), n = k))
          all_distances[[k]] <- distance
        }
        optimal_k <- which.min(all_distances) #which value of k gave closest estimate of forecasted point
        all_optimal_k_values <- optimal_k
      }
    }
  }

  k <- round(mean(all_optimal_k_values)) #what is the average value of k after cross validation

  #############
  #FORECASTING#
  #############
  forecast <- c()
  for (h in 1:forecasting_horizon){
    j <- length(components[[1]])
    all_component_projections <- c()
    for (n in 1:length(components)){
      component <- components[[n]]
      delay_vector <- tail(component, n = m)
      component_projection <- mean(tail(delay_vector, n = k)) #k nearest prediction for component
      all_component_projections <- c(all_component_projections,
                                    component_projection)
      component <- c(component, component_projection) #add prediction to component
      components[[n]] <- component
    }
    x_n_plus_h <- sum(all_component_projections)
    forecast <- c(forecast, x_n_plus_h)
  }

  forecast <- ts(forecast,
                start = time(training_data)[[length(training_data)]] +
                        1/frequency,
                frequency = frequency)
  return(forecast)
}


#' Prophet Forecast
#'
#' @description Applies facebooks prophet algorithm to make forecasts
#' @export
#' @inheritParams mean_forecast

prophet_forecasting <- function(training_data, forecasting_horizon, frequency,
                                apply_limits = TRUE){

  #format data for prophet
  y <- as.vector(training_data)
  times <- as.vector(time(training_data))
  x <- c()
  if (frequency == 12){
    for (i in times){
      string_i <- as.character(i)
      if (nchar(string_i) == 4){
        string_i <- paste0(string_i, '.0')
      }
      year <- strsplit(string_i, split = '.', fixed = TRUE)[[1]][[1]]
      #convert month expressed as decimal to number in the year
      month_number <- strsplit(string_i, split = '.', fixed = TRUE)[[1]][[2]]
      month_number <- as.numeric(paste0('0.', month_number))
      month <- round((month_number*12 + 1), digits = 0) #have to add 1 because january is classed as 0
      month <- as.character(month)
      date <- paste(year, month, '1', sep = '-')
      x <- c(x, date)
    }
  } else{
    stop('Prophet foreacsting only implemented for monthly data')
  }

  df <- as.data.frame(cbind(x, y))
  names(df) <- c('ds', 'y')

  m <- prophet::prophet(df, interval.width = 0.95) #fits prophet model
  if (frequency == 12){
    future <- prophet::make_future_dataframe(m, periods = forecasting_horizon,
                                    freq = 'month') #df for future projections
  } else{
    stop('Prophet foreacsting only implemented for monthly data')
  }
  prediction <- predict(m, future) #do forecast

  residuals <- as.vector(training_data) -
               head(prediction$yhat, nrow(prediction) - forecasting_horizon)

  output <- list(tail(prediction$yhat, forecasting_horizon),
                tail(prediction$yhat_lower, forecasting_horizon),
                tail(prediction$yhat_upper, forecasting_horizon),
                tail(prediction, forecasting_horizon),
                residuals)

  start_date <- time(training_data)[[length(training_data)]] + 1/frequency
  #convert back to time series object
  for (i in 1:3){
    output[[i]] <- ts(output[[i]],
                    start = start_date,
                    frequency = frequency)
  }
  output[[5]] <- ts(output[[5]],
                  start = start_date,
                  frequency = frequency)

  if (apply_limits){
    output <- apply_limits(output)
  }
  return(output)
}

#' Writes the model output to file
#'
#' @description takes a list of outputs from do_forecasts() which summarise all forecasting results and writes them to file
#' @export
#' @param output_list list of forecast summaries produced by do_forecasts()
#' @param output_file name of output file 

write_output <- function(output_list, output_file){

  write.table(rbind(output_list), file = output_file, append = TRUE,
              row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ',')
}

#' Wrtites the headers of an output file
#'
#' @description writes the headers to output file
#' @export
#' @param output_file name of output file 

write_output_file_headers <- function(output_file){


  headers <- list('Organism', 'Location','Ab', 'Forecasting_Horizon',
                  'Ljung-Box_test_P_value', 'Mean_of_Training_Data',
                  'Average_Trend_of_Training_Data','Variance_of_Training_Data')

  for (forecast_name in names(some_forecasting_functions)){
    headers[[length(headers) + 1]] <- paste0(forecast_name, '_ME')
    headers[[length(headers) + 1]] <- paste0(forecast_name, '_RMSE')
    headers[[length(headers) + 1]] <- paste0(forecast_name, '_sMAPE')
    headers[[length(headers) + 1]] <- paste0(forecast_name, '_PI_Width')
    headers[[length(headers) + 1]] <- paste0(forecast_name, '_PI_Accuracy')
    headers[[length(headers) + 1]] <- paste0(forecast_name, '_residuals_mean')
    headers[[length(headers) + 1]] <- paste0(forecast_name,
                                            '_Portmanteau_test_p_value')
  }

  write.table(headers, file = output_file, row.names = FALSE, col.names = FALSE,
    sep = ',', quote = FALSE)
}


#' Perform time series decomposition on monthly data
decompose_data = function(data, s.window=NULL, ...) {
  if (is.null(s.window)) {
    timespan = max(year(data$Date)) - min(year(data$Date))
    if (timespan > 20) {
      s.window = 19
    } else {
      s.window = "periodic"
    }
  }
  
  data_ts = ts(data$logCases, start=1, frequency=12)
  
  decomp = stl(data_ts, s.window, ...)$time.series %>%
    as_tibble() %>%
    mutate(Date = data$Date,
           Year = year(Date),
           Month = month(Date),
           Cases = data$Cases,
           trend.seasonal = exp(trend + seasonal),
           seasonal.remainder = exp(seasonal + remainder),
           trend = exp(trend),
           seasonal = exp(seasonal)) %>%
    pivot_longer(cols = -c(Date, Year, Month))
}

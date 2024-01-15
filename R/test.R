source("R/constants.R")

myfun = function(x) {
  withTimeout({
    Sys.sleep(x)
    
    return(x)
  },
  timeout=2,
  onTimeout = "warning")
}


with_timeout(myfun(1), cpu=1, elapsed=1)

withTimeout(myfun(1.1), timeout=1, onTimeout="warning")


sequence = c(1, 3, 5)

results = lapply(sequence, myfun)

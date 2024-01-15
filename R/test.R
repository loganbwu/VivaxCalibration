source("R/constants.R")

myfun = function(x) {
  Sys.sleep(2)
  
  return(x)
}

myfun(1)

with_timeout(myfun(1), cpu=1, elapsed=1)

withTimeout(myfun(1.1), timeout=1, onTimeout="warning")

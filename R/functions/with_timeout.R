#' Execute function with a maximum time
#' 
#' @param cpu maximum time in seconds
#'
#' https://stackoverflow.com/a/53018594
with_timeout <- function(expr, cpu = Inf, elapsed = Inf){
  if (cpu == Inf & elapsed == Inf) {
    stop("Both cpu and elapsed are Inf. Set at least one.")
  }
  expr <- substitute(expr)
  envir <- parent.frame()
  setTimeLimit(cpu = cpu, elapsed = elapsed, transient = TRUE)
  on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE))
  eval(expr, envir = envir)
}

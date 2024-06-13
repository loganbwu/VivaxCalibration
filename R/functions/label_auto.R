#' Label with abbreviations
label_auto = function(x) {
  max_x = max(x, na.rm=T)
  if (max_x > 5e6) {
    paste(x/1e6, "M")
  } else if (max_x > 5e3) {
    paste(x/1e3, "k")
  } else {
    x
  }
}

#' Label with abbreviations but choose suffix per label
label_auto2 = function(x) {
  case_when(x >= 1e6 ~ paste(x/1e6, "M"),
            x >= 1e3 ~ paste(x/1e3, "k"),
            TRUE ~ as.character(x))
}

#' Label with abbreviations but choose suffix per label
#' Allow fractions
label_auto3 = function(x) {
  case_when(x >= 1e6 ~ paste(x/1e6, "M"),
            x >= 1e3 ~ paste(x/1e3, "k"),
            x < 0.5 ~ paste0("1/", format(round(1/x, 0), scientific=F, nsmall=0, trim=TRUE)),
            TRUE ~ as.character(x))
}

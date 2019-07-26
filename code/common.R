library(pacman)
p_load("shiny")
p_load("knitr")
p_load("kableExtra")
p_load("data.table")
p_load("dplyr")
p_load("ggplot2")
p_load("plotly")
p_load("ggpubr")
p_load("foreach")
p_load("glue")
p_load("WGCNA")
p_load("impute")
p_load("limma")
p_load("itertools")
p_load("parallel")
p_load("doSNOW")
source("../code/plots.R")
source("../rsrc/DistributedLmFit.R", chdir=TRUE)
source("../rsrc/utils.R", chdir=TRUE)
set.seed(1234)

knitr::opts_knit$set(
  root.dir = getwd()
)

knitr::opts_chunk$set(
  context   = "render",
  echo      = FALSE, 
  include   = FALSE, 
  out.width = "100%", 
  fig.width = 10, 
  fig.height= 6)


RECOMPUTE <<- FALSE

touch <- function() {
  RECOMPUTE <<- TRUE
}

cache <- function(condition=RECOMPUTE, foo=foo, fname, verbose=FALSE, ...) {
  # If condition evaluates to TRUE, then 
  # recompute using the call to function
  if (condition | !file.exists(fname)) {
    if(verbose)
      message("Recomputing.")
    res <- foo(...)
    saveRDS(res, fname)
    RECOMPUTE <<- FALSE
  } else {
    if (verbose)
      message("Using cached result.")
    res <- readRDS(fname)
  }
  return (res)
}

#learning logistic regression

## Project
rm(list = ls()) #clear env
graphics.off()  #close plots
cat("\014") #clear console

require(foreign)
require(nnet)
require(ggplot2)
require(reshape2)


ml <- read.dta("http://www.ats.ucla.edu/stat/data/hsbdemo.dta")
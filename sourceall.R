
library("compiler")
library("Rcpp")
library("RcppArmadillo")
library("xtable")
library("here")

str = c("./library/misc/", "./library/stats/", "./library/rcpp")
file.sources <- list.files(str, pattern = "*.R$", full.names = TRUE, ignore.case = TRUE)
sapply(file.sources, source)

file.sources <- list.files(str, pattern = "*.cpp$", full.names = TRUE, ignore.case = TRUE)
sapply(file.sources, Rcpp::sourceCpp)

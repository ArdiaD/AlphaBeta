# This file installs the packages required for the project

str <- c(
  "xts",
  "zoo",
  "xtable",
  "moments",
  "sandwich",
  "Rcpp",
  "RcppEigen",
  "pracma",
  "RcppArmadillo",
  "here"
)

for (i in 1:length(str)) {
  install.packages(str[i])
}
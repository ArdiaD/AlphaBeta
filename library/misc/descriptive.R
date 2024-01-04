
f_stats <- function(Y) {
  tot_T   <- dim(Y)[1]
  tot_n   <- dim(Y)[2]
  tmp     <- apply(!is.na(Y), 2, sum)
  stat_T  <- quantile(tmp, probs = c(0.50, 0.10, 0.90), na.rm = TRUE)
  tmp     <- apply(!is.na(Y), 1, sum)
  stat_n  <- quantile(tmp, probs = c(0.50, 0.10, 0.90), na.rm = TRUE)
  tmp     <- 12 * apply(Y, 2, mean, na.rm = TRUE)
  stat_mu <- quantile(tmp, probs = c(0.50, 0.10, 0.90), na.rm = TRUE)
  tmp     <- sqrt(12) * apply(Y, 2, sd, na.rm = TRUE)
  stat_sd <- quantile(tmp, probs = c(0.50, 0.10, 0.90), na.rm = TRUE)
  tmp     <- sqrt(12) * apply(Y, 2, mean, na.rm = TRUE) / apply(Y, 2, sd, na.rm = TRUE)
  stat_sr <- quantile(tmp, probs = c(0.50, 0.10, 0.90), na.rm = TRUE)
  tmp     <- apply(Y, 2, moments::skewness, na.rm = TRUE)
  stat_sk <- quantile(tmp, probs = c(0.50, 0.10, 0.90), na.rm = TRUE)
  tmp     <- apply(Y, 2, moments::kurtosis, na.rm = TRUE)
  stat_ku <- quantile(tmp, probs = c(0.50, 0.10, 0.90), na.rm = TRUE)
  tmp     <- 1 * apply(Y, 2, function(x) quantile(x, probs = 0.10, na.rm = TRUE))
  qtl_10  <- quantile(tmp, probs = c(0.50, 0.10, 0.90), na.rm = TRUE)
  tmp     <- 1 * apply(Y, 2, function(x) quantile(x, probs = 0.90, na.rm = TRUE))
  qtl_90  <- quantile(tmp, probs = c(0.50, 0.10, 0.90), na.rm = TRUE)

  out  <- c(tot_n, tot_T, stat_T, stat_n, stat_mu, stat_sd, stat_sr, stat_sk, stat_ku, qtl_10, qtl_90)
  str1 <- c("tot_n", "tot_T")
  str2 <- paste(rep(c("T", "n", "mu", "sd", "sr", "sk", "ku", "qtl10", "qtl90"), each = 3),
                rep(c("med", "10", "90"), 9), sep = "_")
  names(out) <- c(str1, str2)
  out
}

f_stats_fact_1 <- function(X) {
  tbl <- c()
  for (i in 1:ncol(X)) {
    Xi  <- X[,i,drop = FALSE]
    tmp <- f_stats(Xi)
    tmp <- tmp[c("T_med", "mu_med", "sd_med", "sr_med", "sk_med", "ku_med", "qtl10_med", "qtl90_med")]
    tbl <- rbind(tbl, round(tmp, 2))
  }
  rownames(tbl) <- colnames(X)
  colnames(tbl) <- c("T", "Av", "Sd", "SR", "Sk", "Ku", "10", "90")
  tbl
}




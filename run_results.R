# This file replicates the main results of Ardia et al. (202x), Is it alpha or beta? Decomposing hedge 
# fund returns when models are misspecified, Journal of Financial Economics. 
# Please look at the README file for more details.

# Load functions (several warning message may appear)
rm(list = ls())
suppressMessages(source("sourceall.R"))
options(digits = 10)

# Session info
sink(here::here("outputs", "session_info.txt"))
Sys.info()
sessionInfo()
sink()

# Load HF data
load(file = here::here("data", "Db_HF_Pseudo.rda")) 
#load(file = here::here("data", "Db_HF.rda")) # This one cannot be shared
names(Db_HF)

# Excess returns
Y <- Db_HF$ret

## Chi1 and chi2 settings
ctr <- list(chi1_max = 15, chi2_max = nrow(Y) / 60)

## DATA STATS / Table 1 -------------------------------------------------------------------------------

strat    <- Db_HF$strat
# strat    <- Db_HF$substrat # !!! Uncomment for subcategories !!!
ustrat   <- sort(unique(strat))
ustrat   <- c(ustrat, "ALL") # Add ALL so that we have global analysis
n_ustrat <- length(ustrat)

tbl1 <- c()
for (i in 1:n_ustrat) {
      cat(ustrat[i], "\n")
      idx  <- f_idx_strat(strat, ustrat[i])
      Yi   <- Y[,idx]
      tmp  <- f_chi(Yi, ctr = ctr)
      idx  <- tmp$sel_Y
      ni   <- sum(idx)
      Yew  <- apply(Yi, 1, mean, na.rm = TRUE) # use full sample
      Yew  <- as.matrix(Yew)
      tmp  <- f_stats(Yew)
      tmp  <- tmp[c("mu_med", "sd_med", "sk_med", "ku_med", "qtl10_med", "qtl90_med")]
      tbl1 <- rbind(tbl1, c(ni, round(tmp, 2)))
}
rownames(tbl1) <- ustrat
colnames(tbl1) <- c("N", "Mean(Ann.)", "Sd(Ann.)", "Sk", "Ku", "10(Ann.)", "90(Ann.)")
tbl1 <- tbl1[,c(2,3,4,5,6,7)]
tbl1
write.table(tbl1, file = here::here("outputs", "Table1.txt"))

# FACTOR MODELS -------------------------------------------------------------------------------

# Load factor data
load(file = here::here("data", "Db_Factors_Pseudo.rda"))
#load(file = here::here("data", "Db_Factors.rda")) # This one cannot be shared

# Build factor models
names(Factors)
X_CAPM <- as.matrix(Factors[,c("FF_EMKT")])
X_FF4  <- as.matrix(Factors[,c("FF_EMKT", "FF_HML", "FF_SMB", "FF_MOM")])
X_FF5  <- as.matrix(Factors[,c("FF_EMKT", "FF_HML", "FF_SMB", "FF_CMA", "FF_RMW")])
X_FH   <- as.matrix(Factors[,c("FF_EMKT", "FF_SMB", "FH_TERM_FH", "FH_DFLT_FH", "FH_PTFSBD", "FH_PTFSCOM", "FH_PTFSFX")])
X_AMP  <- as.matrix(Factors[,c("FF_EMKT", "AQR_VME_VAL", "AQR_VME_MOM")])
X_JKKT <- as.matrix(Factors[,c("FF_EMKT", "FF_HML", "FF_SMB", "FF_MOM", "AQR_TSMOM", "PS_LIQV", "AQR_BAB_USA")])
X_CP   <- as.matrix(Factors[,c("FF_EMKT", "FF_SMB", "PS_LIQV", "AQR_BAB_USA", "BH_VARVIX", "KMPV_GCF", "AQR_TSMOM")])
X_KNS1 <- as.matrix(Factors[,c("FF_EMKT", "KNS_dindrrevlv", "KNS_dindmomrev", "KNS_dindrrev", "KNS_dseason", "KNS_dsue")])
X_KNS2 <- as.matrix(Factors[,c("FF_EMKT", "KNS_pc_d5", "KNS_pc_d1", "KNS_pc_d6", "KNS_pc_d12", "KNS_pc_d2")])

# Gather in list for simplicity
l_NM <- list(
      X_CAPM   = X_CAPM,
      X_FF4    = X_FF4,
      X_FF5    = X_FF5,
      X_FH     = X_FH,
      X_AMP    = X_AMP,
      X_KNS1   = X_KNS1,
      X_KNS2   = X_KNS2,
      X_JKKT   = X_JKKT,
      X_CP     = X_CP
)

## FACTORS / Table 2 -------------------------------------------------------------------------------

nam1 <- c("FF_EMKT", "FF_SMB", "FF_HML", "FF_MOM", "FF_CMA", "FF_RMW", 
          "AQR_VME_VAL", "AQR_VME_MOM",
          "FH_TERM_FH", "FH_DFLT_FH",
          "FH_PTFSBD", "FH_PTFSCOM", "FH_PTFSFX")
nam2 <- c("PS_LIQV", "AQR_BAB_USA",
          "BH_VARVIX",
          "KMPV_GCF", "AQR_TSMOM")
nam  <- c(nam1, nam2)
X <- Factors[,nam]
X[,"BH_VARVIX"] <- X[,"BH_VARVIX"]/10 # Divide by 10 in the table
(tbl1 <- f_stats_fact_1(X))

nam1_ <- c("Market", "Size", "Value", "Momentum", "Investment", "Profitability",
           "Value Everywhere", "Momentum Everywhere",
           "Term", "Default",
           "Bond Straddle", "Commodity Straddle", "Currency Straddle")
nam2_ <- c("Iliquidity", "Betting Against Beta", "Variance", 
           "Carry", "Time-Series Momentum")
nam_  <- c(nam1_, nam2_)
rownames(tbl1) <- nam_
tbl1 <- tbl1[,c("Av", "Sd", "Sk", "Ku", "10", "90")]
tbl1
write.table(tbl1, file = here::here("outputs", "Table2.txt"))

l_nam <- list(nam1, nam2)
n_nam <- length(l_nam)
tbl <- matrix(NA, ncol = n_nam, nrow = n_nam)
for (i in 1:n_nam) {
      for (j in i:n_nam) {
            R <- cor(X[,l_nam[[i]]], X[,l_nam[[j]]])
            tbl[i,j] <- mean(abs(R[lower.tri(R)]))
      }
}

nam <- c("Standard Factors", "Additional Factors")
dimnames(tbl) <- list(nam, nam)
print(round(tbl, 2))

R <- cor(X, use = "pairwise.complete.obs")
R <- R[lower.tri(R)] 
R <- abs(R)
sum(R > 0.5)
length(R)

## ESTIMATION MODELS -------------------------------------------------------------------------------

nrow(Y)
l_M <- l_NM; nam_M <- "NM"

n_M <- length(l_M)
n_Y <- ncol(Y)

m_alpha <- m_R2 <- matrix(data = NA, nrow = n_Y, ncol = n_M)

for (i in 1:n_M) {
      cat("factor model", names(l_M)[i], "\n")
      X   <- l_M[[i]]
      fit <- f_ols(Y = Y, X = X, ctr = ctr)
      print(sum(fit$sel_Y))
      m_alpha[,i] <- fit$alpha_hat
      m_R2[,i]    <- fit$R2_adj
}
colnames(m_alpha) <- colnames(m_R2) <- names(l_M)
rownames(m_alpha) <- rownames(m_R2) <- Db_HF$id

## R2 / Table 4 -------------------------------------------------------------------------------

l_out <- list()
for (i in 1:n_ustrat) {
      idx <- f_idx_strat(strat, ustrat[i])
      tbl <- c()
      for (j in 1:n_M) {
            perf     <- m_alpha[idx,j]
            perf     <- 12 * perf[!is.na(perf)]
            perf_mu  <- mean(perf) # !!! DA later we should have NA
            perf_sd  <- sd(perf)
            perf_neg <- 100 * mean(perf < 0)
            perf_pos <- 100 * mean(perf >= 0)
            perf_q10 <- quantile(perf, probs = 0.10)
            perf_q90 <- quantile(perf, probs = 0.90)
            
            R2      <- 100 * m_R2[idx,j]
            R2      <- R2[!is.na(R2)]
            R2_mean <- mean(R2)
            R2_sd   <- sd(R2)
            R2_q10  <- quantile(R2, probs = 0.1)
            R2_q90  <- quantile(R2, probs = 0.9)
            
            tbl <- rbind(tbl, round(c(perf_mu, perf_sd, 
                                      perf_q10, perf_q90, 
                                      perf_pos, perf_neg,  
                                      rho = NA,
                                      R2_mean, R2_sd, R2_q10, R2_q90),2))
      }
      
      colnames(tbl) <- c("Mean(Ann.)", "Std(Ann.)", 
                         "Q10(Ann.)", "Q90(Ann.)",
                         "Pos", "Neg", 
                         "Rho(CAPM)",
                         "AvR2", "SdR2", "Q10R2", "Q90R2")
      print(ustrat[i])
      rownames(tbl) <- names(l_M)
      l_out[[i]] <- tbl[,c(1,2,6,5,3,4,8,9,10,11)]
}
names(l_out) <- ustrat

tbl <- l_out[[c("ALL")]]
tbl
write.table(tbl, file = here::here("outputs", "Table4r2.txt"))

## AC AND BC / Table 6 -------------------------------------------------------------------------------

f_tbl_ac <- function() {
      tbl_ac <- c()
      tbl_ac_num <- c()
      for (i in 1:length(l_M)) {
            cat("Model ", names(l_M)[i], "\n")
            tmp1 <- f_ctrb(Y = Y[,pos], X_k = l_M[[i]], X_l = NULL, ctr = ctr, type = type, do_diff = FALSE)
            tmp2 <- unlist(tmp1)[c("mean_hat", "vol_hat", "pneg_hat", "ppos_hat", "Q10_hat", "Q90_hat")]
            tmp3 <- unlist(tmp1)[c("mean_se", "vol_se", "pneg_se", "ppos_se", "Q10_se", "Q90_se")]
            
            sc <- c(12, 12, 100, 100, 12, 12)
            tmp2 <- sc * tmp2
            tmp3 <- sc * tmp3
            tbl_ac     <- rbind(tbl_ac, paste0(format(round(tmp2, 2), nsmall = 2), 
                                               " (", format(round(tmp3, 2), nsmall = 2), ")"))
            tbl_ac_num <- rbind(tbl_ac_num, tmp2)
      }
      rownames(tbl_ac) <- names(l_NM)
      colnames(tbl_ac) <- c("Mean", "Vol", "Neg", "Pos", "Q10", "Q90")
      print(tbl_ac)
      sum(pos)
      return(tbl_ac)
}

f_tbl_bc <- function() {
      tbl_bc <- c()
      tbl_bc_num <- c()
      for (i in 1:length(l_M)) {
            cat("Model ", names(l_M)[i], "\n")
            tmp1 <- f_ctrb(Y = Y[,pos], X_k = l_M[[i]], ctr = ctr, type = "beta")
            tmp2 <- unlist(tmp1)[seq(from = 1, to = 12, by = 2)]
            tmp3 <- unlist(tmp1)[seq(from = 2, to = 12, by = 2)]
            
            sc <- c(12, 12, 100, 100, 12, 12)
            tmp2 <- sc * tmp2
            tmp3 <- sc * tmp3
            tbl_bc     <- rbind(tbl_bc, paste0(format(round(tmp2, 2), nsmall = 2), 
                                               " (", format(round(tmp3, 2), nsmall = 2), ")"))
            tbl_bc_num <- rbind(tbl_bc_num, tmp2)
      }
      rownames(tbl_bc) <- names(l_NM)
      colnames(tbl_bc) <- c("Mean", "Vol", "Neg", "Pos", "Q10", "Q90")
      print(tbl_bc)
      sum(pos)
      return(tbl_bc)
}

pos <- rep(TRUE, ncol(Y)); l_M <- l_NM; nam_M <- "NM"
# Uncomment below to generate Table 8 in full
# pos <- f_idx_strat(strat, "Agg Equity 2019"); l_M <- l_NM; nam_M <- "NM" 
# pos <- f_idx_strat(strat, "Agg Macro 2019");  l_M <- l_NM; nam_M <- "NM"
# pos <- f_idx_strat(strat, "Agg Arbitrage 2019");  l_M <- l_NM; nam_M <- "NM"

sum(pos) # Number of funds in category
tmp <- f_chi(Y[,pos], X = l_M[[1]], ctr = ctr); print(sum(tmp$sel_Y)) # Number of funds selected

type <- "alpha"
tbl_ac <- f_tbl_ac()
nam <- c("CAPM", "Carhart", "Five-Factor", "Fung-Hsieh", "AMP", "KNS1", "KNS2", "JKKT", "CP")  
rownames(tbl_ac) <- nam
tbl_ac
write.table(tbl_ac, file = here::here("outputs", "Table6a.txt"))

tbl_bc <- f_tbl_bc()
rownames(tbl_bc) <- nam
tbl_bc
write.table(tbl_bc, file = here::here("outputs", "Table6b.txt"))

## AC DIFFERENCE / Table 5 -----------------------------------------------------------

f_tbl_diff_ac <- function(nam_k, nam_l, do_print = FALSE) {
      tmp3 <- f_ctrb(Y = Y[,pos], X_k = l_M[[nam_k]], X_l = l_M[[nam_l]], 
                     ctr = ctr, type = type, do_diff = TRUE, do_print = do_print)
      tmp4 <- unlist(tmp3)[c("mean_hat", "vol_hat", "pneg_hat", "ppos_hat", "Q10_hat", "Q90_hat")]
      tmp5 <- unlist(tmp3)[c("mean_se", "vol_se", "pneg_se", "ppos_se", "Q10_se", "Q90_se")]
      tmp6 <- abs(tmp4 / tmp5)
      
      sc <- c(12, 12, 100, 100, 12, 12)
      tmp4 <- sc * tmp4
      tmp5 <- sc * tmp5
      
      str <- rep("", length(tmp6))
      str[tmp6 > abs(qnorm(0.05))]  <- "$^{*}$"
      str[tmp6 > abs(qnorm(0.025))] <- "$^{**}$"
      str[tmp6 > abs(qnorm(0.005))] <- "$^{***}$"
      
      tbl1 <- paste0(format(round(tmp4, 2), nsmall = 2), str, 
                     " (", format(round(tmp5, 2), nsmall = 2), ")")
      tbl_print <- print(paste0(paste0(tbl1, collapse = " & "), "\\"))
      tbl_print
      out <- list(tbl = tmp4, tbl_print = tbl_print)
      out
}

# With CAPM as benchmark model
tmp1 <- f_tbl_diff_ac(nam_k = "X_FF4", nam_l = "X_CAPM")
tmp2 <- f_tbl_diff_ac(nam_k = "X_FF5", nam_l = "X_CAPM")
tmp3 <- f_tbl_diff_ac(nam_k = "X_FH", nam_l = "X_CAPM")
tmp4 <- f_tbl_diff_ac(nam_k = "X_AMP", nam_l = "X_CAPM")
tmp5 <- f_tbl_diff_ac(nam_k = "X_KNS1", nam_l = "X_CAPM")
tmp6 <- f_tbl_diff_ac(nam_k = "X_KNS2", nam_l = "X_CAPM")
tmp7 <- f_tbl_diff_ac(nam_k = "X_JKKT", nam_l = "X_CAPM")
tmp8 <- f_tbl_diff_ac(nam_k = "X_CP", nam_l = "X_CAPM")

tbl_diff_ac <- rbind(tmp1$tbl_print, tmp2$tbl_print, tmp3$tbl_print, tmp4$tbl_print, 
                     tmp5$tbl_print, tmp6$tbl_print, tmp7$tbl_print, tmp8$tbl_print)
rownames(tbl_diff_ac) <- nam[-1]
tbl_diff_ac
write.table(tbl_diff_ac, file = here::here("outputs", "Table5.txt"))

## %AC AND BC / Table 4 -------------------------------------------------------------------------------

n <- length(l_M)
tbl <- matrix(NA, nrow = n, ncol = 8)
dimnames(tbl) <- list(names(l_M), c("RETS", "CTR_ALPHA", "CTR_MKT", "CTR_OTHERS", "TOTAL", 
                                    "%CTR_ALPHA", "%CTR_MKT", "%CTR_OTHERS"))
for (i in 1:n) {
      cat(names(l_M)[i], "\n")
      tmp <- f_ctrb_fast(Y[,pos], X_k = l_M[[i]], ctr = ctr)
      tbl[i,1] <- 12 * mean(tmp$mean_ret)
      tbl[i,2] <- 12 * mean(tmp$alpha, na.rm = TRUE)
      tbl[i,3] <- 12 * mean(tmp$bc[,1], na.rm = TRUE)
      tbl[i,4] <- 12 * mean(apply(tmp$bc[,-1], 1, sum), na.rm = TRUE)
      tbl[i,5] <- sum(tbl[i,2:4])
      abstot <- tmp$mean_ret
      abstot[abstot == 0] <- NA
      tbl[i,6] <- 100 * mean(tmp$alpha / abstot, na.rm = TRUE)
      tbl[i,7] <- 100 * mean(tmp$bc[,1] / abstot, na.rm = TRUE)
      tbl[i,8] <- 100 * mean(apply(tmp$bc[,-1], 1, sum) / abstot, na.rm = TRUE)
}

tbl <- round(tbl[,6:8], 2)
rownames(tbl) <- c("CAPM", "Carhart", "Five-Factor", "Fung-Hsieh", 
                   "AMP", "KNS1", "KNS2", "JKKT", "CP") 
tbl
write.table(tbl, file = here::here("outputs", "Table4.txt"))

## INDIVIDUAL BC / Table 7a -------------------------------------------------------------------------------

f_tbl_bc_ind <- function() {
      sc <- c(12, 12, 100, 100, 12, 12)
      tmp2 <- tmp[seq(from = 1, to = 12, by = 2)]
      tmp2 <- matrix(unlist(tmp2), nrow = length(tmp2$mean_hat), byrow = FALSE)
      av_tmp2 <- apply(tmp2[-c(1,nrow(tmp2)),,drop = FALSE], 2, mean)
      tmp2 <- matrix(sc, nrow = nrow(tmp2), ncol = ncol(tmp2), byrow = TRUE) * tmp2
      nam <- names(tmp$mean_hat)
      dimnames(tmp2) <- list(nam, c("Mean", "Vol", "Neg", "Pos", "Q10", "Q90"))
      
      tmp3 <- tmp[seq(from = 2, to = 12, by = 2)]
      tmp3 <- matrix(unlist(tmp3), nrow = length(tmp3$mean_se), byrow = FALSE)
      tmp3 <- matrix(sc, nrow = nrow(tmp3), ncol = ncol(tmp3), byrow = TRUE) * tmp3
      dimnames(tmp3) <- list(nam, c("Mean", "Vol", "Neg", "Pos", "Q10", "Q90"))
      
      tbl <- paste0(format(round(tmp2, 2), nsmall = 2), " (", format(round(tmp3, 2), nsmall = 2), ")")
      tbl <- matrix(tbl, nrow = nrow(tmp2), ncol = ncol(tmp2), byrow = FALSE)
      tbl <- rbind(tbl, format(round(sc * av_tmp2, 2), nsmall = 2))
      dimnames(tbl) <- list(c(nam, "Average"), c("Mean", "Vol", "Neg", "Pos", "Q10", "Q90"))
      
      print(xtable::xtable(tbl))
      sum(pos)
      return(tbl)
}

colnames(l_M[["X_CP"]])
tmp <- f_ctrb_ind_beta(Y = Y[,pos], X = l_M[["X_CP"]], ctr = ctr)
tbl <- f_tbl_bc_ind()
rownames(tbl) <- c("Market", "Size", "Illiquidity", "Betting Against Beta", 
                   "Variance", "Carry", "Time-Series Momentum", 
                   "Total (Others)", "Average (w/o Market)")
tbl
write.table(tbl, file = here::here("outputs", "Table7a.txt"))

## CORRELATIONS BC FOR CP MODEL in Table 7b -------------------------------------------------------------------------------

X_k <- l_M[["X_CP"]]
tmp <- f_ctrb_fast(Y, X_k = X_k, ctr = ctr)
names(tmp)
dim(tmp$bc)

R <- cor(cbind(tmp$alpha, tmp$bc), use = "complete.obs")
dimnames(R) <- list(c("alpha", colnames(X_k)), c("alpha", colnames(X_k)))
tbl <- round(R, 2)
tbl <- tbl[-c(1,8),-c(1,2)]
tbl
write.table(tbl, file = here::here("outputs", "Table7b.txt"))

## DYNAMIC ANALYSIS / Figures 3a-4 -------------------------------------------------------------------------------
## !!! This piece of code is slow !!!

dim(Y)
time <- Db_HF$Date
length(time)

# NB: when using MF, I am also removing the first returns
n_roll <- 205 # 1994-2020
X <- l_M[["X_CP"]] # for hedge funds
# X <- l_M[["X_FF4"]] # for mutual funds
X_b1 <- l_M[["X_CAPM"]] # for benchmark
X_b2 <- l_M[["X_FH"]] # for hedge funds

m_ret <- matrix(data = NA, nrow = ncol(Y[,pos]), ncol = n_roll)
dimnames(m_ret) <- list(1:ncol(Y[,pos]), 1:n_roll)
m_ac <- m_bc_mkt <- m_bc_nmkt <- m_ret
m_ac_b1 <- m_ac_b2 <- m_ret

a_ctb <- array(data = NA, dim = c(ncol(Y[,pos]), ncol(X), n_roll))
dimnames(a_ctb) <- list(1:ncol(Y[,pos]), colnames(X), 1:n_roll)

m_F <- matrix(data = NA, nrow = ncol(X), ncol = n_roll)
dimnames(m_F) <- list(colnames(X), 1:n_roll)
m_beta <- m_F

for (k in 1:n_roll) {
      cat(k, "\n")
      idx <- (1:(120 + (k-1))) # expanding
      Yi <- Y[idx,pos]
      Xi <- X[idx,]
      m_F[,k] <- apply(Xi, 2, mean)
      
      # CP or FF4
      ctri <- list(chi1_max = 15, chi2_max = nrow(Yi) / 60)
      chi  <- f_chi(Y = Yi, X = Xi, ctr = ctri)
      
      m_ret[chi$sel_Y,k] <- apply(Yi[,chi$sel_Y], 2, mean, na.rm = TRUE)
      
      tmp <- f_ctrb_fast(Y = Yi, X_k = Xi, ctr = ctri)
      m_ac[,k]      <- tmp$alpha
      m_beta[,k]    <- apply(tmp$beta, 1, mean, na.rm = TRUE)
      m_bc_mkt[,k]  <- apply(tmp$bc[,1,drop=FALSE], 1, sum)
      m_bc_nmkt[,k] <- apply(tmp$bc[,-1,drop=FALSE], 1, sum)
      a_ctb[,,k]    <- tmp$bc
      
      # Benchmark (CAPM)
      Xi  <- X_b1[idx,,drop=FALSE]
      chi <- f_chi(Y = Yi, X = Xi, ctr = ctri)
      tmp <- f_ctrb_fast(Y = Yi, X_k = Xi, ctr = ctri)
      m_ac_b1[,k] <- tmp$alpha
      
      # Benchmark (FH)
      Xi  <- X_b2[idx,,drop=FALSE]
      chi <- f_chi(Y = Yi, X = Xi, ctr = ctri)
      tmp <- f_ctrb_fast(Y = Yi, X_k = Xi, ctr = ctri)
      m_ac_b2[,k] <- tmp$alpha
}

# Save the results 
save(m_ac, m_beta, m_F, m_bc_mkt, m_bc_nmkt, a_ctb, m_ac_b1, m_ac_b2, 
     file = here::here("outputs", "res_dyn_hf.rda"))

# Load results
load(file = here::here("outputs", "res_dyn_hf.rda"))

# Average
tmp1 <- apply(m_ac, 2, mean, na.rm = TRUE)
tmp2 <- apply(m_bc_mkt, 2, mean, na.rm = TRUE)
tmp3 <- apply(m_bc_nmkt, 2, mean, na.rm = TRUE)
tmp <- 12 * cbind(tmp1, tmp2, tmp3)
colnames(tmp) <- c("Alpha", "Beta Mkt", "Beta Non-Mkt")

ylim <- c(-2, 8)
pdf(file = here::here("outputs", "Figure3a.pdf"), width = 11, height = 7)
par(mfrow = c(1,1))
matplot(tmp, type = "l", lwd = 2, col = c("black", "red", "blue"), 
        lty = c("solid", "dashed", "dotdash"), 
        las = 1, axes = FALSE, ylab = "", ylim = ylim)
axis(side = 2, at = seq(from = -2, to = 8, by = 1), las = 1, cex.axis = 1.3, tck = 1, 
     lty = "dotted", col = "grey")
axis(side = 1, at = seq(from = 1, to = n_roll, by = 12), las = 2,
     labels = seq(from = 2003, to = 2020, by = 1), cex.axis = 1.3)
abline(h = 0)
box()
legend("topright", cex = 1.3, lwd = 2, bg = "white",
       legend = c("Alpha", "Beta Mkt", "Beta Non-Mkt"),
       col = c("black", "red", "blue"), 
       lty = c("solid", "dashed", "dotdash"))
dev.off()

# Benchmarks differences
tmp1_b1 <- apply(m_ac_b1, 2, mean, na.rm = TRUE)
tmp1_d1 <- tmp1_b1 - tmp1 # CAPM - CP
tmp1_b2 <- apply(m_ac_b2, 2, mean, na.rm = TRUE)
tmp1_d2 <- tmp1_b2 - tmp1 # FH - CP
tmp1_b3 <- tmp1_b1 - tmp1_b2
tmp_d   <- 12 * cbind(tmp1_d1, tmp1_d2, tmp1_b3)
colnames(tmp_d) <- c("CAPMvsCP", "FHvsCP", "FHvsCAPM")
tmp_d   <- tmp_d[,c("FHvsCP", "CAPMvsCP", "FHvsCAPM")]

pdf(file = here::here("outputs", "Figure4.pdf"), width = 11, height = 7)
par(mfrow = c(1,1))
matplot(tmp_d, type = "l", lwd = 2, col = c("black", "red", "blue"), 
        lty = c("solid", "dashed", "dotdash"), 
        las = 1, axes = FALSE, ylab = "", ylim = 12 * c(0, 0.4))
axis(side = 2, at = seq(from = 0, to = 5, by = 0.5), las = 1, cex.axis = 1.3, tck = 1, 
     lty = "dotted", col = "grey")
axis(side = 1, at = seq(from = 1, to = n_roll, by = 12), las = 2,
     labels = seq(from = 2003, to = 2020, by = 1), cex.axis = 1.3)
abline(h = 0)
box()
legend("topright", cex = 1.3, lwd = 2, bg = "white",
       legend = c("Fung-Hsieh vs CP", 
                  "CAPM vs CP", 
                  "CAPM vs Fung-Hsieh"),
       col = c("black", "red", "blue"), 
       lty = c("solid", "dashed", "dotdash"))
dev.off()

## FLOWS ANALYSIS / Table 11 -------------------------------------------------------------------------------

AUM <- Db_HF$aum
RET <- Y
dates <- Db_HF$dates
dim(AUM)
dim(RET)

l_dates <- list()
l_dates[[1]] <- list(start = as.Date("1996-01-01"), end = as.Date("2000-12-31"))
l_dates[[2]] <- list(start = as.Date("2001-01-01"), end = as.Date("2005-12-31"))
l_dates[[3]] <- list(start = as.Date("2006-01-01"), end = as.Date("2010-12-31"))
l_dates[[4]] <- list(start = as.Date("2011-01-01"), end = as.Date("2015-12-31"))
l_dates[[5]] <- list(start = as.Date("2016-01-01"), end = as.Date("2020-12-31"))

tmp <- matrix(NA, nrow = 5, ncol = 5)
dimnames(tmp) <- list(c("1996-2000", 
                        "2001-2005", 
                        "2006-2010", 
                        "2011-2015", 
                        "2016-2020"), 
                      c("Low", "2", "3", "4", "High"))
m_nb <- tmp 
m_flow_its <- m_ret_its <- m_ac_its <- tmp
m_ac_pct_its <- m_bc1_pct_its <- m_bc_pct_its <- tmp

store_ret_l <- store_ret_h <- store_ac_l <- store_bc_l <- store_ac_h <- store_bc_h <- NULL
RET_l <- RET_h <- NULL
v_OK_AUM <- rep(NA, 5)

l_RET <- vector("list", 5)

for (i in 1:5) {
      cat("i ", i, "\n")
      pos_its <- dates >= l_dates[[i]]$start & dates <= l_dates[[i]]$end
      pos_its <- which(pos_its)
      
      AUM_its <- AUM[pos_its,]
      RET_its <- RET[pos_its,] 
      X_its   <- l_NM[["X_CP"]][pos_its,]
     
      pos1 <- which(apply(!is.na(RET_its), 2, all))
      pos2 <- which(apply(!is.na(AUM_its), 2, all))
      pos  <- intersect(pos1, pos2)
      v_OK_AUM[i] <- length(pos) / length(pos1)
      
      print(length(pos))
      
      AUM_its   <- AUM_its[,pos]
      RET_its   <- RET_its[,pos]
      FLOWS_its <- (AUM_its[2:60,] / AUM_its[1:59,]) - (1 + RET_its[2:60,] / 100)
      FLOWS_its <- 100 * apply(FLOWS_its, 2, mean)
      
      qtl <- quantile(FLOWS_its, c(0.00, 0.20, 0.40, 0.60, 0.80, 1.00))
      
      for (j in 1:5) {
            cat("j: ", j, "\n")
            pos <- which((FLOWS_its >= qtl[j]) & (FLOWS_its < qtl[j+1]))
            m_nb[i,j] <- length(pos)
            
            tmp <- f_ctrb_fast(RET_its[,pos], X_k = X_its, ctr = ctr)
            m_flow_its[i,j]  <- mean(FLOWS_its[pos])
            
            if (j == 1) {
                  tmp_RET <- matrix(NA, nrow = nrow(RET), ncol = m_nb[i,j])
                  tmp_RET[pos_its,] <- RET_its[,pos]
                  RET_l <- cbind(RET_l, tmp_RET)
            }
            if (j == 5) {
                  tmp_RET <- matrix(NA, nrow = nrow(RET), ncol = m_nb[i,j])
                  tmp_RET[pos_its,] <- RET_its[,pos]
                  RET_h <- cbind(RET_h, tmp_RET)
            }
            
            tmp_RET <- matrix(NA, nrow = nrow(RET), ncol = m_nb[i,j])
            tmp_RET[pos_its,] <- RET_its[,pos]
            l_RET[[j]] <- cbind(l_RET[[j]], tmp_RET)
      }
}

f_tbl <- function() {
      sc <- c(12, 12, 100, 100, 12, 12)
      tmp2 <- tmp[seq(from = 1, to = 12, by = 2)]
      tmp2 <- matrix(unlist(tmp2), nrow = length(tmp2$mean_hat), byrow = FALSE)
      tmp2 <- matrix(sc, nrow = nrow(tmp2), ncol = ncol(tmp2), byrow = TRUE) * tmp2
      nam <- names(tmp$mean_hat)
      dimnames(tmp2) <- list(nam, c("Mean", "Vol", "Neg", "Pos", "Q10", "Q90"))
      
      tmp3 <- tmp[seq(from = 2, to = 12, by = 2)]
      tmp3 <- matrix(unlist(tmp3), nrow = length(tmp3$mean_se), byrow = FALSE)
      tmp3 <- matrix(sc, nrow = nrow(tmp3), ncol = ncol(tmp3), byrow = TRUE) * tmp3
      dimnames(tmp3) <- list(nam, c("Mean", "Vol", "Neg", "Pos", "Q10", "Q90"))
      
      tbl <- paste0(format(round(tmp2, 2), nsmall = 2), " (", format(round(tmp3, 2), nsmall = 2), ")")
      tbl <- matrix(tbl, nrow = nrow(tmp2), ncol = ncol(tmp2), byrow = FALSE)
      dimnames(tbl) <- list(nam, c("Mean", "Vol", "Neg", "Pos", "Q10", "Q90"))
      
      print(xtable::xtable(tbl))
      return(tbl)
}

l_tbl <- list()
for (i in 1:5) {
      cat(i, "\n")
      tmp  <- f_ctrb(Y = l_RET[[i]], X_k = l_NM[["X_CP"]], ctr = ctr, type = "alpha")
      tbl1 <- f_tbl()
      tmp  <- f_ctrb_ind_beta(Y = l_RET[[i]], X = l_NM[["X_CP"]], ctr = ctr)
      tbl2 <- f_tbl()
      tbl2 <- tbl2[c(1,nrow(tbl2)),]
      tbl  <- rbind(tbl1, tbl2)
      rownames(tbl) <- paste(i, c("AC", "BC(Mkt)", "Low BC(Others)"))
      l_tbl[[i]] <- tbl
}
print(l_tbl)

(tbl1 <- t(sapply(l_tbl, function(x) x[1,])))
(tbl2 <- t(sapply(l_tbl, function(x) x[2,])))
(tbl3 <- t(sapply(l_tbl, function(x) x[3,])))

write.table(tbl1, file = here::here("outputs", "Table11a.txt"))
write.table(tbl2, file = here::here("outputs", "Table11b.txt"))
write.table(tbl3, file = here::here("outputs", "Table11c.txt"))

## DISTRIBUTION ALPHA / Figure 2 -------------------------------------------------------------------------------

tmp <- f_ctrb_fast(Y, X_k = l_M[["X_CP"]], ctr = ctr)
ac_cp <- 12 * tmp$alpha; ac_cp <- ac_cp[!is.na(ac_cp)]
bc_cp <- 12 * apply(tmp$bc, 1, sum); bc_cp <- bc_cp[!is.na(bc_cp)]

tmp <- f_ctrb_fast(Y, X_k = l_M[["X_JKKT"]], ctr = ctr)
ac_jkkt <- 12 * tmp$alpha; ac_jkkt <- ac_jkkt[!is.na(ac_jkkt)]
bc_jkkt <- 12 * apply(tmp$bc, 1, sum); bc_jkkt <- bc_jkkt[!is.na(bc_jkkt)]

tmp <- f_ctrb_fast(Y, X_k = l_M[["X_FH"]], ctr = ctr)
ac_fh <- 12 * tmp$alpha; ac_fh <- ac_fh[!is.na(ac_fh)]
bc_fh <- 12 * apply(tmp$bc, 1, sum); bc_fh <- bc_fh[!is.na(bc_fh)]

tmp <- f_ctrb_fast(Y, X_k = l_M[["X_CAPM"]], ctr = ctr)
ac_capm <- 12 * tmp$alpha; ac_capm <- ac_capm[!is.na(ac_capm)]
bc_capm <- 12 * apply(tmp$bc, 1, sum); bc_capm <- bc_capm[!is.na(bc_capm)]

xl <- -20 ; xh <- 20; ylim <- c(0, 0.20) # yearly figure
pdf(file = here::here("outputs", "Figure2a.pdf"), width = 11, height = 7)
tmp <- f_np_pdf(ac_capm, h = NULL, n_mesh = 500, lower = xl, upper = xh, type = 1)
plot(tmp$x, tmp$pdf, xlim = c(xl, xh), type = "l", lwd = 2, ylim = ylim, las = 1, 
     ylab = "", xlab = "", axes = FALSE, lty = "dashed", col = "blue")
mtext("Alpha Component", side = 1, line = 2.2, cex = 1.3)
axis(side = 1, at = seq(from = xl, to = xh, by = 5), cex.axis = 1.3)
axis(side = 2, las = 1, cex.axis = 1.3)
tmp <- f_np_pdf(ac_jkkt, h = NULL, n_mesh = 500, lower = xl, upper = xh, type = 1)
lines(tmp$x, tmp$pdf, col = "red", lwd = 2, lty = "dotdash")
tmp <- f_np_pdf(ac_cp, h = NULL, n_mesh = 500, lower = xl, upper = xh, type = 1)
lines(tmp$x, tmp$pdf, col = "black", lwd = 2, lty = "solid")
abline(v = 0); abline(h = 0); box()
legend("topright", cex = 1.3, lwd = 2, 
       legend = c("CAPM", "JKKT Model", "CP Model"),
       col = c("blue", "red", "black"), 
       lty = c("dashed", "dotdash", "solid"))
dev.off()

pdf(file = here::here("outputs", "Figure2b.pdf"), width = 11, height = 7)
tmp <- f_np_pdf(bc_capm, h = NULL, n_mesh = 500, lower = xl, upper = xh, type = 1)
plot(tmp$x, tmp$pdf, xlim = c(xl, xh), type = "l", lwd = 2, ylim = ylim, las = 1, 
     ylab = "", xlab = "", axes = FALSE, lty = "dashed", col = "blue")
mtext("Beta Component", side = 1, line = 2.2, cex = 1.3)
axis(side = 1, at = seq(from = xl, to = xh, by = 5), cex.axis = 1.3)
axis(side = 2, las = 1, cex.axis = 1.3)
tmp <- f_np_pdf(bc_jkkt, h = NULL, n_mesh = 500, lower = xl, upper = xh, type = 1)
lines(tmp$x, tmp$pdf, col = "red", lwd = 2, lty = "dotdash")
tmp <- f_np_pdf(bc_cp, h = NULL, n_mesh = 500, lower = xl, upper = xh, type = 1)
lines(tmp$x, tmp$pdf, col = "black", lwd = 2, lty = "solid")
abline(v = 0); abline(h = 0); box()
legend("topright", cex = 1.3, lwd = 2, 
       legend = c("CAPM", "JKKT Model", "CP Model"),
       col = c("blue", "red", "black"), 
       lty = c("dashed", "dotdash", "solid"))
dev.off()

## CHARACTERISTICS / Table 10 -------------------------------------------------------------------------------

tmp <- f_chi(Y, X = l_NM[["X_CP"]], ctr = ctr)
print(sum(tmp$sel_Y))  
pos <- which(tmp$sel_Y)

RET <- Y[,pos]
dim(RET)

mf   <- Db_HF$mf[pos]
pf   <- Db_HF$pf[pos]
hwm  <- Db_HF$hwm[pos]
hr   <- Db_HF$hr[pos]
np   <- Db_HF$np[pos] # in months
lp   <- Db_HF$lp[pos] # in years

l_RET <- list()
l_RET[[1]] <- RET[, which(mf >= 0 & mf < 2.0)]
l_RET[[2]] <- RET[, which(mf >= 2.0)]
names(l_RET) <- c("<2.0", ">=2.0")

# !!! Uncomment below to generate other outputs in Table !!!
# l_RET <- list()
# l_RET[[1]] <- RET[, which(pf >= 0 & pf < 20)]
# l_RET[[2]] <- RET[, which(pf >= 20)]
# names(l_RET) <- c("<20", ">=20")
# 
# l_RET <- list()
# l_RET[[1]] <- RET[, which(hwm == TRUE)]
# l_RET[[2]] <- RET[, which(hwm == FALSE)]
# names(l_RET) <- c("Yes", "No")
# 
# l_RET <- list()
# l_RET[[1]] <- RET[, which(lp == 0)]
# l_RET[[2]] <- RET[, which(lp > 0)]
# names(l_RET) <- c("=0y", ">0y")
# 
# l_RET <- list()
# l_RET[[1]] <- RET[, which(np == 0)]
# l_RET[[2]] <- RET[, which(np > 0)]
# names(l_RET) <- c("=0m", ">0m")

n_RET <- length(l_RET)
l_tbl <- vector('list', n_RET)
for (i in 1:n_RET) {
      cat(i, "\n")
      tmp  <- f_ctrb(Y = l_RET[[i]], X_k = l_NM[["X_CP"]], ctr = ctr, type = "alpha")
      l_tbl[[i]] <- f_tbl()
}
print(l_tbl)

tbl1 <- t(sapply(l_tbl, function(x) x[1,]))
rownames(tbl1) <- names(l_RET)
tbl1
write.table(tbl1, file = here::here("outputs", "Table10.txt"))

## TOY EXAMPLE / Figure 1 -------------------------------------------------------------------------------

mu_a   <- 0.0
sig2_a <- 1.5^2
lam    <- 7.5
mu_b   <- 0.3
sig2_b <- 0.4^2  

mu_b1 <- 0.3 #mu_b
mu_b2 <- 0.1 #mu_b / 3
mu_b3 <- 0.1 #mu_b / 3

f_ac_mu_0 <- function()
      mu_a + lam * (mu_b1 + mu_b2 + mu_b3)

f_ac_sig2_0 <- function()
      sig2_a + 3 * lam^2 * sig2_b

f_ac_pdf_0 <- function(x)
      dnorm(x, mean = f_ac_mu_0(), sd = sqrt(f_ac_sig2_0()))

f_ac_mu_1 <- function()
      mu_a + lam * mu_b3

f_ac_sig2_1 <- function()
      sig2_a + lam^2 * sig2_b 

f_ac_pdf_1 <- function(x)
      dnorm(x, mean = f_ac_mu_1(), sd = sqrt(f_ac_sig2_1()))

pdf(file = here::here("outputs", "Figure1a.pdf"), width = 11, height = 7)
x <- seq(from = -20, to = 20, length.out = 1000)
plot(x, f_ac_pdf_0(x), type = "l", las = 1, ylim = c(0, 0.15), 
     col = "blue", lwd = 2, lty = "dashed", ylab = "", 
     xlab = "", axes = FALSE)
lines(x, f_ac_pdf_1(x), type = "l", col = "black", lwd = 2, lty = "solid")
mtext("Alpha Component", side = 1, line = 2.2, cex = 1.3)
axis(side = 1, at = seq(from = -20, to = 20, by = 5), cex.axis = 1.3)
axis(side = 2, las = 1, cex.axis = 1.3)
legend("topright", legend = c("CAPM", "HF Model"), cex = 1.3,
       lwd = c(2, 2), col = c("blue", "black"), lty = c("dashed","solid"))
abline(v = 0); box()
dev.off()

f_bc_mu_0 <- function()
      mu_b * lam

f_bc_sig2_0 <- function()
      sig2_b * lam^2

f_bc_pdf_0 <- function(x)
      dnorm(x, mean = f_bc_mu_0(), sd = sqrt(f_bc_sig2_0()))

f_bc_mu_1 <- function()
      (mu_b + mu_b1 + mu_b2) * lam

f_bc_sig2_1 <- function()
      3 * lam^2 * sig2_b

f_bc_pdf_1 <- function(x)
      dnorm(x, mean = f_bc_mu_1(), sd = sqrt(f_bc_sig2_1()))

pdf(file = here::here("outputs", "Figure1b.pdf"), width = 11, height = 7)
x <- seq(from = -20, to = 20, length.out = 1000)
plot(x, f_bc_pdf_0(x), type = "l", las = 1, ylim = c(0, 0.15), 
     col = "blue", lwd = 2, lty = "dashed", ylab = "", 
     xlab = "", axes = FALSE)
lines(x, f_bc_pdf_1(x), type = "l", col = "black", lwd = 2, lty = "solid")
mtext("Beta Component", side = 1, line = 2.2, cex = 1.3)
axis(side = 1, at = seq(from = -20, to = 20, by = 5), cex.axis = 1.3)
axis(side = 2, las = 1, cex.axis = 1.3)
legend("topright", legend = c("CAPM", "HF Model"), cex = 1.3,
       lwd = c(2, 2), col = c("blue", "black"), lty = c("dashed","solid"))
abline(v = 0); box()
dev.off()



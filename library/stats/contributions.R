
f_winztrunc <- function(x, type, p) {
      
      x_ <- x
      
      if (type == "winz1") {
            ax   <- abs(x)
            aqtl <- as.numeric(quantile(ax, prob = 1 - p), na.rm = TRUE)
            x_[ax >= aqtl]  <- aqtl
            x_[ax <= -aqtl] <- -aqtl
      }
      
      if (type == "winz2") { # winzorization two-sided
            qtl <- as.numeric(quantile(x, prob = c(p/2, 1 - p/2)), na.rm = TRUE)
            x_[x_ <= qtl[1]] <- qtl[1]
            x_[x_ >= qtl[2]] <- qtl[2]
      }
      
      if (type == "trim1") {
            ax   <- abs(x)
            aqtl <- as.numeric(quantile(ax, prob = 1 - p), na.rm = TRUE)
            x_[ax >= aqtl]  <- 0
            x_[ax <= -aqtl] <- 0
      }
      
      if (type == "trim2") { # trimming two-sided
            qtl <- as.numeric(quantile(x, prob = c(p/2, 1 - p/2)), na.rm = TRUE)
            x_[x_ <= qtl[1]] <- 0
            x_[x_ >= qtl[2]] <- 0
      }
      
      return(x_)
}

f_ctrb_fast <- function(Y, X_k, X_l = NULL, ctr = list(), do_diff = FALSE) {
      
      Y     <- f_Y(Y)
      n_Y   <- ncol(Y)
      
      n_X_k <- ncol(X_k)
      fit_k <- f_ols(Y, X_k, ctr = ctr)
      sel_Y <- fit_k$sel_Y
      
      if (do_diff) {
            fit_l <- f_ols(Y, X_l, ctr = ctr)
            sel_Y <- sel_Y & fit_l$sel_Y
      }
      
      a_hat_k   <- fit_k$alpha_hat
      sig_hat_k <- apply(fit_k$resid, 2, sd, na.rm = TRUE)
      ir_hat_k  <- a_hat_k / sig_hat_k
      
      if (do_diff) {
            a_hat_l   <- fit_l$alpha_hat
            sig_hat_l <- apply(fit_l$resid, 2, sd, na.rm = TRUE)
            ir_hat_l  <- a_hat_l / sig_hat_l
            
            a_hat_k  <- a_hat_k - a_hat_k
            ir_hat_k <- ir_hat_k - ir_hat_l
      }
      
      av_ret   <- apply(Y[,sel_Y,drop = FALSE], 2, mean, na.rm = TRUE)
      mean_ret <- mean(av_ret)
      
      EF <- matrix(data = NA, nrow = n_Y, ncol = n_X_k)
      for (i in 1:n_Y) {
            X_Ri <- X_k
            X_Ri[is.na(Y[,i]),] <- NA
            EF[i,] <- apply(X_Ri, 2, mean, na.rm = TRUE)
      }
      
      b_hat_k <- fit_k$beta_hat 
      bc_k    <- t(b_hat_k) * EF
      
      mean_ac <- mean(a_hat_k, na.rm = TRUE)
      names(mean_ac) <- "alpha"
      
      mean_bc <- apply(bc_k, 2, mean, na.rm = TRUE)
      names(mean_bc) <- colnames(X_k)
      
      out <- list(av_ret   = av_ret,   # time-series average return of selected funds
                  mean_ret = mean_ret, # overall average return of selected funds  
                  mean_ac  = mean_ac,  # overall alpha contribution
                  mean_bc  = mean_bc,  # overall beta contribution
                  alpha    = a_hat_k,  # alpha contribution for all funds
                  ir       = ir_hat_k, # ir of universe
                  beta     = b_hat_k,  # beta_k exposure for all funds
                  bc       = bc_k)     # beta_k contribution for all funds
      return(out)
}

f_ctrb <- function(Y, 
                   X_k, 
                   X_l      = NULL, 
                   ctr      = list(),
                   type     = c("alpha", "ir", "beta"),
                   do_diff  = FALSE,
                   do_print = FALSE) {
      
      f_qtl_se_bahadur <- function(level) {
            Q_k <- as.numeric(quantile(ctrb_k, level))
            Q_l <- 0
            
            if (do_diff) {
                  Q_l <- as.numeric(quantile(ctrb_l, level))
            }
            Q <- Q_k - Q_l 
            fQ_k <- fit_np_k$pdf[which.min(abs(fit_np_k$x - Q_k))]
            tmp_Q_k <- matrix(h_k^(-1) * dnorm((ctrb_k - Q_k) / h_k) * fQ_k^(-1), nrow = nchi, ncol = T) * tmp_k
            
            tmp_Q_l <- 0
            if (do_diff) {
                  fQ_l <- fit_np_l$pdf[which.min(abs(fit_np_l$x - Q_l))]
                  tmp_Q_l <- matrix(h_l^(-1) * dnorm((ctrb_l - Q_l) / h_l) * fQ_l^(-1), nrow = nchi, ncol = T) * tmp_l
            }
            
            tmp <- tmp_Q_k - tmp_Q_l
            tmp <- f_tcrossprod_eigen(tmp)
            tmp <- f_winztrunc(tmp, type = type_wt, p = p_wt)
            V_Q <- sum(tmp) / (nchi^2 * T) 
            out <- list(Q = Q, V_Q = V_Q)
            return(out)
      }
      
      f_qtl_se_inverse <- function(level) {
            
            Q_k <- as.numeric(quantile(ctrb_k, level))
            Q_l <- 0
            
            if (do_diff) {
                  Q_l <- as.numeric(quantile(ctrb_l, level))
            }
            Q <- Q_k - Q_l 
            
            qtl_b     <- c(0.45, 0.55) # current setup, 10%
            iqrn      <- diff(qnorm(qtl_b))
            tmp_P_Q_k <- matrix(-h_k^(-1) * dnorm((ctrb_k - Q_k) / h_k), nrow = nchi, ncol = T) * tmp_k
            
            tmp_P_Q_l <- 0
            if (do_diff) {
                  tmp_P_Q_l <- matrix(-h_l^(-1) * dnorm((ctrb_l - Q_l) / h_l), nrow = nchi, ncol = T) * tmp_l
            }
            
            tmp       <- tmp_P_Q_k - tmp_P_Q_l
            tmp       <- f_tcrossprod_eigen(tmp)
            tmp       <- f_winztrunc(tmp, type = type_wt, p = p_wt)
            V_P_Q     <- sum(tmp) / (nchi^2 * T) 
            qtl_P_Q   <- qnorm(qtl_b, mean = level, sd = sqrt(V_P_Q))
            qtl_P_Q   <- pmax(pmin(qtl_P_Q, 1.0), 0.0)
            
            if (!do_diff) {
                  pos1 <- min(max(sum(cdf_np_k <= qtl_P_Q[2]), 1), length(cdf_np_k)) 
                  pos2 <- min(max(sum(cdf_np_k <= qtl_P_Q[1]), 1), length(cdf_np_k))
                  V_Q  <- ((fit_np_k$x[pos1] - fit_np_k$x[pos2]) / iqrn)^2
            }
            if (do_diff) {
                  pos1 <- min(max(sum(cdf_np_d <= qtl_P_Q[2]), 1), length(cdf_np_d)) 
                  pos2 <- min(max(sum(cdf_np_d <= qtl_P_Q[1]), 1), length(cdf_np_d))
            }
            V_Q <- ((fit_np_k$x[pos1] - fit_np_k$x[pos2]) / iqrn)^2
            
            out <- list(Q = Q, V_Q = V_Q)
            return(out)
      }
      
      ################################################################################################
      
      type_wt  <- "winz2"
      p_wt     <- 0.01
      type     <- type[1]
      
      ac_k <- ac_l <- 0
      ir_k <- ir_l <- 0
      bc_k <- bc_l <- 0
      
      ################################################################################################
      Y     <- f_Y(Y)
      n_X_k <- ncol(X_k)
      Xf_k  <- f_Xf(T = nrow(Y), X = X_k)
      
      fit_k <- f_ols(Y, X_k, ctr = ctr)
      sel_X <- fit_k$sel_X
      sel_Y <- fit_k$sel_Y
      
      if (do_diff) {
            n_X_l <- ncol(X_l)
            Xf_l  <- f_Xf(T = nrow(Y), X = X_l)
            
            fit_l <- f_ols(Y, X_l, ctr = ctr)
            sel_X <- sel_X & fit_l$sel_X
            sel_Y <- sel_Y & fit_l$sel_Y
      }
      
      Y_sel <- Y[sel_X,sel_Y,drop = FALSE]
      Xf_sel_k <- Xf_k[sel_X,,drop = FALSE]
      if (do_diff) {
            Xf_sel_l <- Xf_l[sel_X,,drop = FALSE]
      }
      
      n     <- ncol(Y)
      nchi  <- ncol(Y_sel)
      T     <- nrow(Y_sel)
      Ti    <- apply(!is.na(Y_sel), 2, sum)
      tau_k <- tau_l <- T / Ti
      
      # misspecified case
      T_se <- ceiling(mean(Ti))
      # well-specified case
      #T_se <- nchi
      
      ac_k <- fit_k$alpha_hat[sel_Y]
      b_k  <- fit_k$beta_hat[,sel_Y,drop = FALSE]
      e_k  <- fit_k$resid[sel_X,sel_Y]
      e_k[is.na(e_k)] <- 0
      
      sig_k <- apply(fit_k$resid[sel_X,sel_Y], 2, sd, na.rm = TRUE)
      ir_k  <- ac_k / sig_k
      
      EF_k <- matrix(data = NA, nrow = nchi, ncol = n_X_k)
      for (i in 1:nchi) {
            X_Ri <- X_k
            X_Ri[is.na(Y_sel[,i]),] <- NA
            EF_k[i,] <- apply(X_Ri, 2, mean, na.rm = TRUE)
      }
      
      # beta contribution for whole set of factors
      bc_k <- t(b_k) * EF_k
      bc_k <- apply(bc_k, 1, sum)
      
      # computation of E_1' * Q_x^-1 * eps_i,t * x_t
      E1invQi_k <- f_E1invQi_rcpp(Y_sel, Xf_sel_k) # nchi x K
      tmp_k     <- E1invQi_k %*% t(Xf_sel_k)
      tmp_k     <- tmp_k * t(e_k)
      tmp_k     <- tmp_k * matrix(tau_k, nrow = nchi, ncol = T)
      
      if (type == "ir") {
            tmp1  <- matrix(sig_k, nrow = nchi, ncol = T)
            tmp2  <- matrix(ir_k / (2 * sig_k^3), nrow = nchi, ncol = T) * (t(e_k^2) - tmp1^2)
            tmp_k <- (1 / tmp1) * tmp_k - tmp2
      }
      if (type == "beta") {
            tmp <- t(Y_sel)
            tmp[is.na(tmp)] <- 0
            tmp_k <- tmp - tmp_k
      }
      
      tmp_l <- 0
      if (do_diff) {
            ac_l <- fit_l$alpha_hat[sel_Y]
            b_l  <- fit_l$beta_hat[,sel_Y,drop = FALSE]
            e_l  <- fit_l$resid[sel_X,sel_Y]
            e_l[is.na(e_l)] <- 0
            
            if (type == "ir") {
                  # in case of information ratio
                  sig_l <- apply(fit_l$resid[sel_X,sel_Y], 2, sd, na.rm = TRUE)
                  ir_l  <- ac_l / sig_l
            }
            if (type == "beta") {
                  EF_l <- matrix(data = NA, nrow = nchi, ncol = n_X_l)
                  for (i in 1:nchi) {
                        X_Ri <- X_l
                        X_Ri[is.na(Y_sel[,i]),] <- NA
                        EF_l[i,] <- apply(X_Ri, 2, mean, na.rm = TRUE)
                  }
                  
                  # recompute with full observations for numerical stability
                  EF_l <- matrix(apply(X_l, 2, mean), nrow = nchi, ncol = n_X_l, byrow = TRUE)
                  
                  # beta contribution for whole set of factors
                  bc_l <- t(b_l) * EF_l
                  bc_l <- apply(bc_l, 1, sum)
            }
            
            # computation of E_1' * Q_x^-1 * eps_i,t * x_t
            E1invQi_l <- f_E1invQi_rcpp(Y_sel, Xf_sel_l) # nchi x K
            tmp_l     <- E1invQi_l %*% t(Xf_sel_l)
            tmp_l     <- tmp_l * t(e_l)
            tmp_l     <- tmp_l * matrix(tau_l, nrow = nchi, ncol = T)
            
            if (type == "ir") {
                  tmp1  <- matrix(sig_l, nrow = nchi, ncol = T)
                  tmp2  <- matrix(ir_l / (2 * sig_l^3), nrow = nchi, ncol = T) * (t(e_l^2) - tmp1^2)
                  tmp_l <- (1 / tmp1) * tmp_l - tmp2
            }
            if (type == "beta") {
                  tmp <- t(Y_sel)
                  tmp[is.na(tmp)] <- 0
                  tmp_l <- tmp - tmp_l
            }
      }
      
      if (type == "alpha") {
            ctrb_k <- ac_k
            ctrb_l <- ac_l
      }
      
      if (type == "ir") {
            ctrb_k <- ir_k
            ctrb_l <- ir_l
      }
      if (type == "beta") {
            ctrb_k <- bc_k
            ctrb_l <- bc_l
      }
      
      ################################################################################################
      # M1
      # for mean, g1(a) = 1
      M1_k <- mean(ctrb_k)
      M1_l <- 0
      if (do_diff) {
            M1_l <- mean(ctrb_l)
      }
      M1   <- M1_k - M1_l
      tmp  <- tmp_k - tmp_l
      tmp  <- f_tcrossprod_eigen(tmp)
      tmp  <- f_winztrunc(tmp, type = type_wt, p = p_wt)
      V_M1 <- sum(tmp) / (nchi^2 * T) 
      
      ################################################################################################
      # M2
      # for std, g1(a) = D2M2 * 2 * a_i + D1M2 
      M2_k     <- sqrt(mean(ctrb_k^2) - M1_k^2)
      D2M2_k   <- (2 * M2_k)^(-1)
      D1M2_k   <- -M1_k / M2_k
      tmp_M2_k <- matrix((D2M2_k * 2 * ctrb_k + D1M2_k), nrow = nchi, ncol = T) * tmp_k
      M2_l     <- 0 
      tmp_M2_l <- 0
      if (do_diff) {
            M2_l     <- sqrt(mean(ctrb_l^2) - M1_l^2)
            D2M2_l   <- (2 * M2_l)^(-1)
            D1M2_l   <- -M1_l / M2_l
            tmp_M2_l <- matrix((D2M2_l * 2 * ctrb_l + D1M2_l), nrow = nchi, ncol = T) * tmp_l
      }
      
      M2   <- M2_k - M2_l
      tmp  <- tmp_M2_k - tmp_M2_l
      tmp  <- f_tcrossprod_eigen(tmp)
      tmp  <- f_winztrunc(tmp, type = type_wt, p = p_wt)
      V_M2 <- sum(tmp) / (nchi^2 * T) 
      
      ################################################################################################
      # KERNEL
      h_k <- 1.06 * (length(ctrb_k)^(-1/5)) * sd(ctrb_k)
      fit_np_k <- f_np_pdf(ctrb_k, h = h_k, n_mesh = 1000, lower = min(ctrb_k), upper = max(ctrb_k), type = 1)
      cdf_np_k <- cumsum(fit_np_k$pdf) / sum(fit_np_k$pdf)
      if (do_diff) {
            h_l <- 1.06 * (length(ctrb_l)^(-1/5)) * sd(ctrb_l)
            fit_np_l <- f_np_pdf(ctrb_l, h = h_l, n_mesh = 1000, lower = min(ctrb_l), upper = max(ctrb_l), type = 1)
            cdf_np_l <- cumsum(fit_np_l$pdf) / sum(fit_np_l$pdf)
            
            ctrb_d <- ctrb_k - ctrb_l
            h_d  <- 1.06 * (length(ctrb_d)^(-1/5)) * sd(ctrb_d)
            fit_np_d <- f_np_pdf(ctrb_d, h = h_d, n_mesh = 1000, lower = min(ctrb_d), upper = max(ctrb_d), type = 1)
            cdf_np_d <- cumsum(fit_np_d$pdf) / sum(fit_np_d$pdf)
      }
      
      ################################################################################################
      # PROPORTIONS
      Pneg_k  <- mean(ctrb_k < 0.0)
      tmp_P_k <- matrix(-h_k^(-1) * dnorm((ctrb_k - 0.0) / h_k), nrow = nchi, ncol = T) * tmp_k
      
      Pneg_l  <- 0
      tmp_P_l <- 0
      if (do_diff) {
            Pneg_l  <- mean(ctrb_l < 0.0)
            tmp_P_l <- matrix(-h_l^(-1) * dnorm((ctrb_l - 0.0) / h_l), nrow = nchi, ncol = T) * tmp_l
      }
      
      Pneg   <- Pneg_k - Pneg_l
      tmp    <- tmp_P_k - tmp_P_l
      tmp    <- f_tcrossprod_eigen(tmp)
      tmp    <- f_winztrunc(tmp, type = type_wt, p = p_wt)
      V_Pneg <- sum(tmp) / (nchi^2 * T) 
      
      Ppos_k <- mean(ctrb_k >= 0.0)
      Ppos_l <- 0
      if (do_diff) {
            Ppos_l <- mean(ctrb_l >= 0.0)
      }
      
      Ppos   <- Ppos_k - Ppos_l
      V_Ppos <- V_Pneg
      
      ################################################################################################
      # QUANTILES
      Q10 <- f_qtl_se_bahadur(level = 0.10)
      Q90 <- f_qtl_se_bahadur(level = 0.90)
      
      if (do_print) {
            cat("n:    ", n, "\n")
            cat("nchi: ", nchi, "\n")
            cat("T:    ", T, "\n")
            cat("Tse:  ", T_se, "\n")
      }
      
      out <- list(
            mean_hat = M1,    mean_se = sqrt(V_M1) / sqrt(T_se), 
            vol_hat  = M2,    vol_se  = sqrt(V_M2) / sqrt(T_se), 
            pneg_hat = Pneg,  pneg_se = sqrt(V_Pneg) / sqrt(T_se), 
            ppos_hat = Ppos,  ppos_se = sqrt(V_Ppos) / sqrt(T_se), 
            Q10_hat  = Q10$Q, Q10_se  = sqrt(Q10$V_Q) / sqrt(T_se), 
            Q90_hat  = Q90$Q, Q90_se  = sqrt(Q90$V_Q) / sqrt(T_se)
      )
      
      pos <- which(unlist(lapply(out, length)) == 0)
      if (length(pos) > 0) {
            out[[pos]] <- NA
      }
      
      return(out)
}

f_ctrb_ind_beta <- function(Y, X, ctr = list()) {
      
      type_wt  <- "winz2"
      p_wt     <- 0.01
      type_qtl <- 1
      
      f_qtl_se <- function(level, type_qtl, tmp_P_Q_tot, tmp_Q_tot) {
            
            Q_hat <- as.numeric(quantile(bc_j, level))
            tmp_Q <- tmp_P_Q <- NULL
            
            if (type_qtl == 1) {
                  fQ  <- fit_np$pdf[which.min(abs(fit_np$x - Q_hat))]
                  tmp_Q <- matrix(h^(-1) * dnorm((bc_j_ - Q_hat) / h) * fQ^(-1), nrow = nchi, ncol = T) * tmp_aij
                  
                  if (j <= n_X) {
                        tmp <- f_tcrossprod_eigen(tmp_Q)
                        if (j > 1) {
                              tmp_Q_tot <- tmp_Q_tot + tmp_Q
                        }
                  } else {
                        tmp <- f_tcrossprod_eigen(tmp_Q_tot)
                  }
                  tmp <- f_winztrunc(tmp, type = type_wt, p = p_wt)
                  V_Q_hat <- sum(tmp) / (nchi^2 * T) 
            }

            if (type_qtl == 2) {
                  # asymptotic variance for se of proportion 
                  tmp_P_Q <- matrix(-h^(-1) * dnorm((bc_j_ - Q_hat) / h), nrow = nchi, ncol = T) * tmp_aij
                  if (j <= n_X) {
                        tmp <- f_tcrossprod_eigen(tmp_P_Q)
                  } else {
                        tmp <- f_tcrossprod_eigen(tmp_P_Q_tot)
                  }
                  tmp <- f_winztrunc(tmp, type = type_wt, p = p_wt)
                  V_P_Q_hat <- sum(tmp) / (nchi^2 * T)
                  print(V_P_Q_hat)
                  if (V_P_Q_hat < 0) {
                        cat("negative variance")
                        h <- 1.06 * (length(bc_j_)^(-1/5)) * sd_bc_j
                        if (j <= n_X) {
                              tmp_P_Q <- matrix(-h^(-1) * dnorm((bc_j_ - Q_hat) / h), nrow = nchi, ncol = T) * tmp_aij
                              tmp <- f_tcrossprod_eigen(tmp_P_Q)
                        } else {
                              tmp <- f_tcrossprod_eigen(tmp_P_Q_tot)
                        } 
                        tmp <- f_winztrunc(tmp, type = type_wt, p = p_wt)
                        V_P_Q_hat <- sum(tmp) / (nchi^2 * T)
                  }
                  tmp_P_Q_tot <- tmp_P_Q_tot + tmp_P_Q
                  qtl_b     <- c(0.45, 0.55) # CI interval
                  qtl_P_Q   <- qnorm(qtl_b, mean = level, sd = sqrt(V_P_Q_hat))
                  qtl_P_Q   <- pmax(pmin(qtl_P_Q, 1.0), 0.0)
                  # invert confidence interval
                  pos1      <- min(max(sum(cdf_np <= qtl_P_Q[2]), 1), length(cdf_np)) # find position in CDF
                  pos2      <- min(max(sum(cdf_np <= qtl_P_Q[1]), 1), length(cdf_np)) # find position in CDF
                  # determine the variance of a normal that yield this confidence interval
                  iqrn      <- diff(qnorm(qtl_b))
                  V_Q_hat   <- ((fit_np$x[pos1] - fit_np$x[pos2]) / iqrn)^2
            }
            
            out <- list(Q_hat       = Q_hat, 
                        V_Q_hat     = V_Q_hat, 
                        tmp_Q       = tmp_Q, 
                        tmp_Q_tot   = tmp_Q_tot,
                        tmp_P_Q     = tmp_P_Q, 
                        tmp_P_Q_tot = tmp_P_Q_tot)
            return(out)
      }
      
      Y   <- f_Y(Y)
      n_Y <- ncol(Y)
      n_T <- nrow(Y)
      n_X <- ncol(X)
      nam <- colnames(X)
      
      # b_i,j^k
      fit      <- f_ols(Y = Y, X = X, ctr = ctr)
      beta_i_j <- fit$beta_hat
      eta_i_j  <- fit$resid
      sel_X    <- fit$sel_X
      sel_Y    <- fit$sel_Y
      
      EF_j <- matrix(data = NA, nrow = n_Y, ncol = n_X)
      for (i in 1:n_Y) {
            X_Ri <- X
            X_Ri[is.na(Y[,i]),] <- NA
            EF_j[i,] <- apply(X_Ri, 2, mean, na.rm = TRUE)
      }
      
      # beta contribution for whole set of factors
      bc_i_j <- t(beta_i_j) * EF_j
      
      Xf       <- f_Xf(T = n_Y, X = X)
      Y_sel    <- Y[sel_X,sel_Y,drop = FALSE]
      Xf_sel   <- Xf[sel_X,,drop = FALSE]
      
      nchi <- ncol(Y_sel)
      n    <- ncol(Y)
      T    <- nrow(Y_sel)
      Ti   <- apply(!is.na(Y_sel), 2, sum)
      tau  <- T / Ti
      
      T_se <- ceiling(mean(Ti))
      #T_se <- nchi # well-specified case
      
      eta_i_j[is.na(eta_i_j)] <- 0
      
      tmpNA <- rep(NA, n_X + 1)
      names(tmpNA) <- c(nam, "TOTAL")
      mean_hat <- vol_hat <- pneg_hat <- ppos_hat <- Q10_hat <- Q90_hat <- tmpNA
      mean_se  <- vol_se  <- pneg_se  <- ppos_se  <- Q10_se  <- Q90_se  <- tmpNA
      
      bc_tot <- tmp_M1_tot <- tmp_M2_tot <- tmp_P_tot <- 0 
      tmp_P_Q_tot10 <- tmp_P_Q_tot90 <- tmp_Q_tot10 <- tmp_Q_tot90 <- 0
      
      for (j in 1:(n_X + 1)) {
            if (j <= n_X) {
                  cat("Contribution Factor ", nam[j], "\n")
                  bc_j    <- bc_i_j[sel_Y,j]
                  EjinvQi <- f_EjinvQi_rcpp(Y_sel, Xf_sel, j + 1)
                  tmp_1   <- EjinvQi %*% t(Xf_sel)
                  tmp_1   <- tmp_1 * t(eta_i_j[,sel_Y])
                  tmp_1   <- matrix(EF_j[sel_Y,j], nrow = nchi, ncol = T) * tmp_1
                  tmp_2   <- matrix(Xf_sel[,j + 1], nrow = nchi, ncol = T, byrow = TRUE)
                  tmp_2   <- matrix(beta_i_j[j,sel_Y], nrow = nchi, ncol = T) * tmp_2
                  tmp_aij <- tmp_1 + tmp_2
                  if (j > 1) {
                        bc_tot <- bc_tot + bc_j # total contribution
                  }
            } else {
                  cat("Total Contribution\n")
                  bc_j <- bc_tot  
            }
            
            ################################################################################################
            cat("mean\n")
            # for mean, g1(a) = 1
            M1_j <- mean(bc_j)
            tmp_M1 <- tmp_aij
            if (j <= n_X) {
                  tmp <- f_tcrossprod_eigen(tmp_M1)
                  if (j > 1) {
                        tmp_M1_tot <- tmp_M1_tot + tmp_M1
                  }
            } else {
                  tmp <- f_tcrossprod_eigen(tmp_M1_tot)
            }
            tmp    <- f_winztrunc(tmp, type = type_wt, p = p_wt)
            V_M1_j <- sum(tmp) / (nchi^2 * T)
            
            ################################################################################################
            cat("std\n")
            # for std, g1(a) = D2M2 * 2 * a_i + D1M2 
            M2_j   <- sqrt(mean(bc_j^2) - M1_j^2)
            D2M2   <- (2 * M2_j)^(-1)
            D1M2   <- -M1_j / M2_j
            tmp_M2 <- matrix((D2M2 * 2 * bc_j + D1M2), nrow = nchi, ncol = T) * tmp_aij
            if (j <= n_X) {
                  tmp <- f_tcrossprod_eigen(tmp_M2)
                  if (j > 1) {
                        tmp_M2_tot <- tmp_M2_tot + tmp_M2
                  }
            } else {
                  tmp <- f_tcrossprod_eigen(tmp_M2_tot)
            }
            tmp    <- f_winztrunc(tmp, type = type_wt, p = p_wt)
            V_M2_j <- sum(tmp) / (nchi^2 * T)
            
            ################################################################################################
            bc_j_   <- bc_j
            sd_bc_j <- sd(bc_j_)
            h  <- 1.06 * (length(bc_j_)^(-1/5)) * sd_bc_j
            LB <- min(bc_j_) #- 0.5 * diff(range(bc_j_))
            UB <- max(bc_j_) #+ 0.5 * diff(range(bc_j_))
            fit_np  <- f_np_pdf(bc_j_, h = h, n_mesh = 1000, lower = LB, upper = UB, type = 1)
            cdf_np  <- cumsum(fit_np$pdf) / sum(fit_np$pdf)
            
            ################################################################################################
            cat("prop\n")
            Pneg_j   <- mean(bc_j < 0.0)
            Ppos_j   <- mean(bc_j >= 0.0)
            tmp_dst  <- -h^(-1) * dnorm((bc_j_ - 0.0) / h)
            tmp_P    <- matrix(tmp_dst, nrow = nchi, ncol = T) * tmp_aij
            if (j <= n_X) {
                  tmp <- f_tcrossprod_eigen(tmp_P)
                  if (j > 1) {
                        tmp_P_tot <- tmp_P_tot + tmp_P
                  }
            } else {
                  tmp <- f_tcrossprod_eigen(tmp_P_tot)
            }
            tmp      <- f_winztrunc(tmp, type = type_wt, p = p_wt)
            V_Pneg_j <- sum(tmp) / (nchi^2 * T)
            V_Ppos_j <- V_Pneg_j
            #print(V_Ppos_j)
            
            ################################################################################################
            cat("qtl\n")
            Q10_j <- f_qtl_se(level       = 0.10, 
                              type_qtl    = type_qtl, 
                              tmp_P_Q_tot = tmp_P_Q_tot10, 
                              tmp_Q_tot   = tmp_Q_tot10)
            
            Q90_j <- f_qtl_se(level       = 0.90, 
                              type_qtl    = type_qtl, 
                              tmp_P_Q_tot = tmp_P_Q_tot90,
                              tmp_Q_tot   = tmp_Q_tot90)
            
            tmp_P_Q_tot10 <- Q10_j$tmp_P_Q_tot
            tmp_P_Q_tot90 <- Q90_j$tmp_P_Q_tot
            tmp_Q_tot10   <- Q10_j$tmp_Q_tot
            tmp_Q_tot90   <- Q90_j$tmp_Q_tot
            
            mean_hat[j] <- M1_j
            vol_hat[j]  <- M2_j
            pneg_hat[j] <- Pneg_j
            ppos_hat[j] <- Ppos_j
            Q10_hat[j]  <- Q10_j$Q_hat
            Q90_hat[j]  <- Q90_j$Q_hat
            
            mean_se[j]  <- sqrt(V_M1_j) / sqrt(T_se) 
            vol_se[j]   <- sqrt(V_M2_j) / sqrt(T_se)
            pneg_se[j]  <- sqrt(V_Pneg_j) / sqrt(T_se)
            ppos_se[j]  <- sqrt(V_Ppos_j) / sqrt(T_se)
            Q10_se[j]   <- sqrt(Q10_j$V_Q_hat) / sqrt(T_se)
            Q90_se[j]   <- sqrt(Q90_j$V_Q_hat) / sqrt(T_se)
      }
      
      out <- list(
            mean_hat = mean_hat, mean_se = mean_se,
            vol_hat  = vol_hat,  vol_se  = vol_se,
            pneg_hat = pneg_hat, pneg_se = pneg_se,
            ppos_hat = ppos_hat, ppos_se = ppos_se,
            Q10_hat  = Q10_hat,  Q10_se  = Q10_se,
            Q90_hat  = Q90_hat,  Q90_se  = Q90_se
      )
      
      pos <- which(unlist(lapply(out, length)) == 0)
      if (length(pos) > 0) {
            out[[pos]] <- NA
      }
      
      return(out)
}



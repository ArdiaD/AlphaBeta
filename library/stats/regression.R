
## f_Y and f_Xf --------------------------------------------------------------------------

f_Y <- function(Y) {
  Y <- as.matrix(Y)
  dimnames(Y) <- NULL
  Y
}

f_Xf <- function(T, X = NULL) {
  if (is.null(X)) {
    Xf <- matrix(rep(1, T), nrow = T, ncol = 1)
  } else {
    X <- as.matrix(X)
    dimnames(X) <- NULL
    Xf <- cbind(1, X)
  }
  Xf
}

## f_chi --------------------------------------------------------------------------------

f_chi_dft <- function(ctr) {
  ctr_dft <- list(chi1_max = Inf, chi2_max = Inf, chi3_min = 0)
  ctr_dft[names(ctr)] <- ctr
  ctr_dft
}

# Selection function
f_chi <- function(Y, X = NULL, ctr = list()) {
  
  if (is.null(ctr)) {
    out <- list(sel_X = rep(TRUE, nrow(Y)),
                sel_Y = rep(TRUE, ncol(Y)))
    return(out)
  }
  
  ctr       <- f_chi_dft(ctr)
  Y         <- f_Y(Y)
  Xf        <- f_Xf(X, T = nrow(Y))
  sel_X     <- apply(!is.na(Xf), 1, all)
  tmp       <- f_chi_rcpp(Y[sel_X,,drop=FALSE], Xf[sel_X,,drop=FALSE], ctr)
  n_obs     <- as.vector(tmp$n_obs)
  CN        <- as.vector(tmp$CN)
  chi0_isok <- as.logical(tmp$chi0_isok)
  chi1_isok <- as.logical(tmp$chi1_isok) # CN
  chi2_isok <- as.logical(tmp$chi2_isok) # T/Ti
  chi3_isok <- as.logical(tmp$chi3_isok) # nobs min
  
  sel_Y <- (chi0_isok & chi1_isok & chi2_isok & chi3_isok)
  
  out <- list(
    sel_X     = sel_X,
    sel_Y     = sel_Y,
    n_obs     = n_obs,
    CN        = CN,
    chi0_isok = chi0_isok,
    chi1_isok = chi1_isok,
    chi2_isok = chi2_isok,
    chi3_isok = chi3_isok
  )
  
  return(out)
}

## f_ols ----------------------------------------------------------------------

f_ols <- function(Y, X = NULL, ctr = list()) {
  
  Y  <- f_Y(Y)
  T  <- nrow(Y)
  Xf <- f_Xf(T, X)
  n  <- ncol(Y)
  k1 <- ncol(Xf)
  k  <- k1 - 1
  
  n_obs     <- rep(NA, n)
  alpha_hat <- alpha_se <- alpha_tstat <- alpha_pval <- rep(NA, n)
  R2        <- R2_adj   <- sig_hat     <- sig_res    <- rep(NA, n)
  beta_hat  <- beta_se  <- beta_tstat  <- beta_pval  <- NULL
  resid     <- matrix(data = NA, nrow = T, ncol = n)
  if (k >= 1) {
    beta_hat <- beta_se <- beta_tstat <- beta_pval <- matrix(data = NA, nrow = k, ncol = n)
  }
  
  # first check the filter
  tmp <- f_chi(Y, X, ctr = ctr)
  sel_X <- tmp$sel_X
  sel_Y <- tmp$sel_Y
  fit <- f_ols_rcpp(Y[sel_X,sel_Y,drop=FALSE], Xf[sel_X,,drop=FALSE])
  #browser()
  alpha_hat[sel_Y] <- as.vector(fit$gamma[1,])
  alpha_se[sel_Y]  <- as.vector(fit$sigXX[1,])
  alpha_tstat      <- alpha_hat / alpha_se
  if (k >= 1) {
    beta_hat[,sel_Y] <- matrix(fit$gamma[2:(k + 1),], nrow = k, ncol = sum(sel_Y), byrow = FALSE)
    beta_se[,sel_Y]  <- matrix(fit$sigXX[2:(k + 1),], nrow = k, ncol = sum(sel_Y), byrow = FALSE)
    beta_tstat       <- beta_hat / beta_se
  }
  
  n_obs[sel_Y]       <- as.vector(fit$n_obs)
  sig_hat[sel_Y]     <- as.vector(fit$sig_hat)
  sig_res[sel_Y]     <- as.vector(fit$sig_res)
  R2[sel_Y]          <- as.vector(fit$R2)
  R2_adj[sel_Y]      <- as.vector(fit$R2_adj)
  resid[sel_X,sel_Y] <- matrix(fit$resid, nrow = sum(sel_X), ncol = sum(sel_Y), byrow = FALSE) 
  
  out <- list(
    n_obs       = n_obs,
    alpha_hat   = alpha_hat, 
    alpha_se    = alpha_se,
    alpha_tstat = alpha_tstat,
    beta_hat    = beta_hat,
    beta_se     = beta_se,
    beta_tstat  = beta_tstat,
    sig_hat     = sig_hat,
    sig_res     = sig_res,
    R2          = R2,
    R2_adj      = R2_adj,
    resid       = resid,
    sel_X       = sel_X,
    sel_Y       = sel_Y
  )
  
  return(out)
}
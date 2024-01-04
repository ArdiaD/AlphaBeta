
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat f_scaleXf_rcpp_old_and_wrong(const arma::mat & Xf) {
      
      arma::mat id(Xf.n_rows, 1, arma::fill::ones);
      arma::mat m;
      arma::mat s;
      arma::mat Xs;
      
      m = arma::mean(Xf, 0);
      s = arma::stddev(Xf, 0);
      m(0,0) = 0;
      s(0,0) = 1;
      
      Xs = (Xf - arma::kron(id, m)) / arma::kron(id, s);
      //Xs = Xf / arma::kron(id, s);
      
      return Xs;
}

// [[Rcpp::export]]
arma::mat f_scaleXf_rcpp(const arma::mat & Xf) {

      arma::mat id(Xf.n_rows, 1, arma::fill::ones);
      arma::mat m;
      arma::mat s;
      arma::mat Xs;

      arma::rowvec norms = arma::sqrt(arma::sum(arma::square(Xf), 0));

      norms.replace(0, 1); // Make sure norms are not zero to avoid division by zero

      Xs = Xf / arma::kron(id, norms);

      return Xs;
}

// [[Rcpp::export]]
double f_sumcrossprod_rcpp(const arma::mat & X) {
  return arma::as_scalar(accu(X.t() * X));
}

// [[Rcpp::export]]
double f_sumcrossprod2_rcpp(const arma::mat & X) {
  
  int n = X.n_cols;
  int T = X.n_rows;
  double sumcross = 0;
  
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int t = 0; t < T; ++t) {
        sumcross = sumcross + X(t,i) * X(t,j);
      }
    }
  }
  
  return sumcross;
}

// [[Rcpp::export]]
arma::mat f_crossprod_rcpp(const arma::mat & X) {
  
  int n = X.n_cols;
  int T = X.n_rows;
  double sumt = 0;
  arma::mat XXt(n, n, arma::fill::zeros);
  
  for (int i1 = 0; i1 < n; ++i1) {
    for (int i2 = 0; i2 < n; ++i2) {
      sumt = 0;
      for (int t = 0; t < T; ++t) {
        sumt = sumt + X(t,i1) * X(t,i2);
      }
      XXt(i1,i2) = sumt;
      XXt(i2,i1) = sumt;
    }
  }
  
  return XXt;
}

// [[Rcpp::export]]
arma::mat f_E1invQi_rcpp(const arma::mat & Y, const arma::mat & Xf) {
  
  int n  = Y.n_cols;
  int k1 = Xf.n_cols;
  
  arma::mat Xi;
  arma::mat XXi;
  arma::mat Qi;
  arma::mat invQi;
  arma::uvec pos;
  int Ti;
  
  arma::mat E1invQi(n, k1, arma::fill::ones);
  
  for (int i = 0; i < n; ++i) {
    pos    = find_finite(Y.col(i));
    Xi     = Xf.rows(pos);
    Ti     = Xi.n_rows;
    XXi    = (Xi.t() * Xi);
    Qi     = XXi / Ti;
    invQi  = inv_sympd(Qi);
    E1invQi.row(i) = invQi.row(0);
  }
  
  return E1invQi;
}

// [[Rcpp::export]]
arma::mat f_EjinvQi_rcpp(const arma::mat & Y, const arma::mat & Xf, const int j) {
  
  int n  = Y.n_cols;
  int k1 = Xf.n_cols;
  
  arma::mat Xi;
  arma::mat XXi;
  arma::mat Qi;
  arma::mat invQi;
  arma::uvec pos;
  int Ti;
  
  arma::mat EjinvQi(n, k1, arma::fill::ones);
  
  for (int i = 0; i < n; ++i) {
    pos    = find_finite(Y.col(i));
    Xi     = Xf.rows(pos);
    Ti     = Xi.n_rows;
    XXi    = (Xi.t() * Xi);
    Qi     = XXi / Ti;
    invQi  = inv_sympd(Qi);
    EjinvQi.row(i) = invQi.row(j-1);
  }
  
  return EjinvQi;
}

// [[Rcpp::export]]
List f_chi_rcpp(const arma::mat & Y, const arma::mat & Xf, List ctr) {
  
  int n  = Y.n_cols;
  int T  = Y.n_rows;
  int k1 = Xf.n_cols;
  
  double chi1_max = ctr["chi1_max"];
  double chi2_max = ctr["chi2_max"];
  double chi3_min = ctr["chi3_min"];
  
  arma::colvec yi;
  arma::mat Xi;
  arma::mat XXi;
  arma::mat Qi;
  arma::uvec pos;
  arma::vec tmp;
  int Ti;
  
  arma::vec n_obs(n, arma::fill::zeros);
  arma::vec CN(n, arma::fill::zeros);
  CN.fill(1000);
  
  arma::vec chi0_isok(n, arma::fill::zeros);
  arma::vec chi1_isok(n, arma::fill::zeros);
  arma::vec chi2_isok(n, arma::fill::zeros);
  arma::vec chi3_isok(n, arma::fill::zeros);
  chi0_isok.fill(false);
  chi1_isok.fill(false);
  chi2_isok.fill(false);
  chi3_isok.fill(false);
  
  for (int i = 0; i < n; ++i) {
    yi  = Y.col(i);
    pos = find_finite(yi);
    Ti  = pos.size();
    n_obs(i) = Ti;
    
    //CN(i) = 1000;
    if (Ti > 0) {
      // chi0 - regression feasible
      if (Ti >= k1) {
        chi0_isok(i) = true;
        
        // condition number
        Xi  = Xf.rows(pos);
        Xi  = f_scaleXf_rcpp(Xi);
        XXi = (Xi.t() * Xi);
        Qi  = XXi / Ti;
        tmp = arma::eig_sym(Qi);
        if (all(tmp > 0)) {
          CN(i) = sqrt(tmp(k1-1) / tmp(0));
        }
      }
      // chi1 - condition number
      if (CN(i) <= chi1_max) {
        chi1_isok(i) = true;
      }
      
      // chi2 - enough data
      if (((double)T / (double)Ti) <= chi2_max) {
        chi2_isok(i) = true;
      }
    }
    // chi3 - min obs
    if (Ti >= chi3_min) {
      chi3_isok(i) = true;
    }
  }
  
  return List::create(Named("n_obs")     = n_obs, 
                      Named("CN")        = CN,
                      Named("chi0_isok") = chi0_isok,
                      Named("chi1_isok") = chi1_isok,
                      Named("chi2_isok") = chi2_isok,
                      Named("chi3_isok") = chi3_isok);
}

// [[Rcpp::export]]
List f_ols_rcpp(const arma::mat & Y, const arma::mat & Xf) {
  
  int n  = Y.n_cols;
  int T  = Y.n_rows;
  int k1 = Xf.n_cols;
  
  arma::colvec yi;
  arma::mat Xi;
  arma::mat XXi;
  arma::mat Qi;
  arma::mat invXXi;
  arma::uvec pos;
  double sig2;
  double rssi;
  int Ti;
  arma::colvec gammai;
  arma::colvec resi;
  arma::mat tmp(T, 1, arma::fill::zeros);
  
  arma::vec n_obs(n, arma::fill::zeros);
  arma::mat gamma(k1, n, arma::fill::zeros); 
  arma::vec sig_hat(n, arma::fill::zeros);
  arma::vec sig_res(n, arma::fill::zeros);
  arma::mat sigXX(k1, n, arma::fill::zeros); 
  arma::vec R2(n, arma::fill::zeros);
  arma::vec R2_adj(n, arma::fill::zeros);
  arma::mat resid(T, n, arma::fill::zeros);
  
  for (int i = 0; i < n; ++i) {
    yi       = Y.col(i);
    pos      = find_finite(yi);
    yi       = yi(pos);
    Ti       = yi.size();
    n_obs(i) = Ti;
    Xi       = Xf.rows(pos);
    XXi      = (Xi.t() * Xi);
    invXXi   = inv_sympd(XXi);
    gammai   = invXXi * Xi.t() * yi;
    resi     = yi - Xi * gammai;
    rssi     = arma::as_scalar(arma::sum(arma::pow(resi,2)));
    sig2     = rssi / (Ti - k1);
    sigXX.col(i) = arma::sqrt(sig2 * arma::diagvec(invXXi));
    
    gamma.col(i) = gammai;
    sig_hat(i)   = std::sqrt(sig2);
    sig_res(i)   = std::sqrt(rssi/ Ti);
    R2(i)        = (1 - rssi / arma::sum( pow(yi - arma::mean(yi), 2) ) );
    R2_adj(i)    = (1 - ((1 - R2(i)) * (Ti - 1)) / (Ti - k1));
    tmp.fill(NA_REAL);
    tmp(pos)     = resi;
    resid.col(i) = tmp;
  }
  
  return List::create(Named("n_obs")   = n_obs, 
                      Named("gamma")   = gamma, 
                      Named("sig_hat") = sig_hat,
                      Named("sig_res") = sig_res,
                      Named("R2")      = R2,
                      Named("R2_adj")  = R2_adj,
                      Named("sigXX")   = sigXX,
                      Named("resid")   = resid);
}


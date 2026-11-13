#include <iostream>
#include <armadillo>
#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

struct MstepResult{
  arma::vec beta;
  arma::vec lambda;
  arma::vec cumulative_lambda;
  arma::vec exp_beta_x;
};

struct EMResult{
  arma::vec beta;
  arma::vec lambda;
  arma::vec cumulative_lambda;
};

static inline arma::vec ComputeExpBetaX(const arma::vec& beta, const arma::mat& X){
  return arma::exp(X * beta);
}

// Time is in descending order; return the index of the first entry <= val, or m if none exists (out-of-range marker)
static inline arma::uword first_leq_desc(const arma::vec& Time, double val){
  arma::uword m = Time.n_elem;
  arma::uword lo = 0, hi = (m==0)?0:(m-1), ans = m;
  if(m==0) return m;
  while(lo <= hi){
    arma::uword mid = (lo + hi) >> 1;
    if(Time(mid) <= val){
      ans = mid;
      if(mid == 0) break;
      hi = mid - 1;
    }else{
      lo = mid + 1;
    }
  }
  return ans;
}

// compute the tail sum from k to the end (k in [0, m]; 0 when k == m)
static inline arma::vec tail_cumsum(const arma::vec& lambda){
  arma::uword m = lambda.n_elem;
  arma::vec tail(m + 1, arma::fill::zeros);
  for(arma::sword k = (arma::sword)m - 1; k >= 0; --k){
    tail((arma::uword)k) = tail((arma::uword)k + 1) + lambda((arma::uword)k);
  }
  return tail; // length m+1, tail[k] = sum_{j=k}^{m-1} lambda[j]; tail[m] = 0
}

static inline arma::vec CumulativeLambda(const arma::vec& lambda){
  // original logic: accumulate from the end toward the front
  arma::uword m = lambda.n_elem;
  arma::vec cumulative_lambda(m, arma::fill::zeros);
  double acc = 0.0;
  for(arma::sword i = (arma::sword)m - 1; i >= 0; --i){
    acc += lambda((arma::uword)i);
    cumulative_lambda((arma::uword)i) = acc;
  }
  return cumulative_lambda;
}

static inline double AbsoluteDifferenceSafe(double prev, double next, double delta) {
  double mn = std::abs(prev);
  if(std::abs(next) < mn) mn = std::abs(next);
  if (mn <= delta) return std::abs(prev - next);
  return std::abs(prev - next) / mn;
}

static inline double VectorAbsoluteDifferenceSafe(const arma::vec& prev, const arma::vec& next, double delta) {
  double mx = 0.0;
  for (arma::uword i = 0; i < prev.n_elem; ++i) {
    double d = AbsoluteDifferenceSafe(prev(i), next(i), delta);
    if (d > mx) mx = d;
  }
  return mx;
}

// ====== Likelihood: use binary search + tail sum instead of arma::find ======
static double likelihood_i_fast(const arma::vec& beta,
                                const arma::vec& lambda,
                                const arma::mat& X,
                                const arma::mat& simpdat,
                                const arma::vec& Time,
                                arma::uword i,
                                const arma::vec& exp_beta_x,   // passed in for reuse
                                const arma::vec& tailLam){     // length m+1
  int id_i = static_cast<int>(simpdat(i,0));
  double L_i = simpdat(i,1);
  double R_i = simpdat(i,2);
  double e_bx = exp_beta_x(id_i - 1);
  
  arma::uword m = Time.n_elem;
  arma::uword idxL = first_leq_desc(Time, L_i);
  arma::uword idxR = first_leq_desc(Time, R_i);
  
  if(std::isinf(R_i)){
    double sumL = (idxL < m) ? tailLam(idxL) : 0.0;
    return - e_bx * sumL;
  }else{
    double sumL = (idxL < m) ? tailLam(idxL) : 0.0;
    double sumR = (idxR < m) ? tailLam(idxR) : 0.0;
    double a = std::exp(- e_bx * sumL);
    double b = std::exp(- e_bx * sumR);
    double diff = a - b;
    if(diff < 1e-300) diff = 1e-300;
    return std::log(diff);
  }
}

static double likelihood_beta_fast(const arma::vec& beta,
                                   const arma::vec& lambda,
                                   const arma::mat& X,
                                   const arma::mat& simpdat,
                                   const arma::vec& Time){
  arma::uword n = X.n_rows;
  arma::vec exp_beta_x = ComputeExpBetaX(beta, X);
  arma::vec tailLam = tail_cumsum(lambda);
  
  double pl = 0.0;
  for(arma::uword i = 0; i < n; ++i){
    pl += likelihood_i_fast(beta, lambda, X, simpdat, Time, i, exp_beta_x, tailLam);
  }
  return pl;
}

// ====== profile_likelihood: update only lambda, avoiding building the full w ======
static arma::vec profile_likelihood_fast(const arma::vec& beta,
                                         const arma::vec& lambda_initial,
                                         const arma::mat& X,
                                         const arma::mat& simpdat,
                                         double eps,
                                         const arma::vec& Time){
  arma::uword n = X.n_rows;
  arma::uword m = Time.n_elem;
  
  arma::vec lambda_hat = lambda_initial;
  arma::vec lambda_old = lambda_initial;
  arma::vec exp_beta_x = ComputeExpBetaX(beta, X);
  
  double error = 100.0;
  int Iter = 0;
  const int MaxIter = 1000;
  const double delta = 0.01;
  
  while (error > eps && Iter < MaxIter){
    arma::vec tailLam = tail_cumsum(lambda_old);
    
    // numerator: colsum_w(k); denominator: risk-set sum exp under the new beta
    arma::vec colsum_w(m, arma::fill::zeros);
    
    // first compute (k_start, k_end) and c_i for each i and accumulate into colsum_w
    for(arma::uword i = 0; i < n; ++i){
      int id_i = static_cast<int>(simpdat(i,0));
      double L_i = simpdat(i,1);
      double R_i = simpdat(i,2);
      if(std::isinf(R_i)) continue;
      
      arma::uword idxL = first_leq_desc(Time, L_i);
      arma::uword idxR = first_leq_desc(Time, R_i);
      int k_start = (idxR < m) ? (int)idxR : -1;
      int k_end   = (idxL > 0) ? (int)(idxL - 1) : -1;
      if(k_start < 0 || k_end < 0 || k_start > k_end) continue;
      
      double e_bx = exp_beta_x(id_i - 1);
      double SL = (idxL < m ? tailLam(idxL) : 0.0) * e_bx;
      double SR = (idxR < m ? tailLam(idxR) : 0.0) * e_bx;
      double denom = 1.0 - std::exp(SL - SR);
      if(denom < 1e-12) denom = 1e-12;
      double c = e_bx / denom;
      
      for(int k = k_start; k <= k_end; ++k){
        colsum_w((arma::uword)k) += lambda_old((arma::uword)k) * c;
      }
    }
    
    // denominator: risk-set sum exp for each k
    arma::vec denom_lambda(m, arma::fill::zeros);
    double acc = 0.0;
    arma::uword idx_star = 0;
    for(arma::uword k = 0; k < m; ++k){
      double t_k = Time(k);
      while(idx_star < n && simpdat(idx_star, 3) >= t_k){
        int id_j = static_cast<int>(simpdat(idx_star,4));
        acc += exp_beta_x(id_j - 1);
        ++idx_star;
      }
      denom_lambda(k) = std::max(acc, 1e-16);
    }
    
    for(arma::uword k = 0; k < m; ++k){
      lambda_hat(k) = colsum_w(k) / denom_lambda(k);
    }
    
    // convergence criterion uses the cumulative lambda (same as the original code)
    arma::vec cum_old = CumulativeLambda(lambda_old);
    arma::vec cum_hat = CumulativeLambda(lambda_hat);
    error = VectorAbsoluteDifferenceSafe(cum_hat, cum_old, delta);
    
    lambda_old = lambda_hat;
    ++Iter;
  }
  
  if(error <= eps){
    Rcpp::Rcout << "Lambda converged in " << Iter << " iterations, Error = " << error << std::endl;
  } else{
    Rcpp::Rcout << "Lambda did not converge!" << std::endl;
  }
  
  return lambda_hat;
}


// [[Rcpp::export]]
Rcpp::List score_test_improved(arma::vec beta_new, arma::mat X_new, arma::vec lambda_old,
                      arma::mat simpdat, double eps, arma::vec Time, double constant){
  Rcpp::List ret;
  long n = X_new.n_rows;
  long p = X_new.n_cols;
  double h = constant/sqrt(n);
  
  arma::mat Dpl(n, p);
  
  double pl;
  arma::vec pl1(p);
  arma::mat pl2(p,p); pl2.zeros();
  arma::mat Info(p,p);
  
  //first order
  for(long j = 0; j < p; ++j){
    arma::vec beta_j1 = beta_new;
    arma::vec beta_j2 = beta_new;
    beta_j1(j) = beta_j1(j) + h;
    beta_j2(j) = beta_j2(j) - h;
    arma::vec lambda_j1 = profile_likelihood_fast(beta_j1, lambda_old, X_new, simpdat, eps, Time);
    arma::vec lambda_j2 = profile_likelihood_fast(beta_j2, lambda_old, X_new, simpdat, eps, Time);
    arma::vec exp_bx_j1 = ComputeExpBetaX(beta_j1, X_new);
    arma::vec tailLam_j1 = tail_cumsum(lambda_j1);
    arma::vec exp_bx_j2 = ComputeExpBetaX(beta_j2, X_new);
    arma::vec tailLam_j2 = tail_cumsum(lambda_j2);
    for(long i = 0; i < n; ++i){
      Dpl(i,j) = (likelihood_i_fast(beta_j1, lambda_j1, X_new, simpdat, Time, i, exp_bx_j1, tailLam_j1) - 
        likelihood_i_fast(beta_j2, lambda_j2, X_new, simpdat, Time, i, exp_bx_j2, tailLam_j2)) / (2*h);
    }
  }
  
  // second order
  arma::vec maximizing_lambda = profile_likelihood_fast(beta_new, lambda_old, X_new, simpdat, eps, Time);
  pl = likelihood_beta_fast(beta_new, lambda_old, X_new, simpdat, Time);
  
  for(long j = 0; j < p; ++j){
    arma::vec beta_j = beta_new;
    beta_j(j) = beta_j(j) + h;
    arma::vec lambda_j = profile_likelihood_fast(beta_j, lambda_old, X_new, simpdat, eps, Time);
    pl1(j) = likelihood_beta_fast(beta_j, lambda_j, X_new, simpdat, Time);
  }
  
  for(long i = 0; i < p; ++i){
    for(long j = 0; j <= i; ++j){
      arma::vec beta_ij = beta_new;
      beta_ij(i) = beta_ij(i) + h;
      beta_ij(j) = beta_ij(j) + h;
      arma::vec lambda_ij = profile_likelihood_fast(beta_ij, lambda_old, X_new, simpdat, eps, Time);
      pl2(i,j) = likelihood_beta_fast(beta_ij, lambda_ij, X_new, simpdat, Time);
      Info(i,j) = (pl - pl1(i) - pl1(j) + pl2(i,j)) / (h * h);
      Info(j,i) = Info(i,j);
    }
  }
  
  double Stat = arma::sum(Dpl.col(p-1))/sqrt(n);
  arma::mat Info1 = arma::trans(Dpl) * Dpl/n;
  arma::mat Info2 = -Info/n;
  
  ret["Stat"] = Stat;
  ret["Dpl"] = Dpl;
  ret["Info1"] = Info1;
  ret["Info2"] = Info2;
  
  return ret;
}


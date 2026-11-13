#include <iostream>
#include <armadillo>
#include <RcppArmadillo.h>
#include <cmath>
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

// ====== Core: M-step without explicitly storing w ======
static MstepResult Mstep_memory_friendly(const arma::vec& beta,
                                         const arma::mat& X,
                                         const arma::mat& simpdat,
                                         const arma::vec& Time,
                                         const arma::vec& lambda){
  
  const arma::uword p = beta.n_elem;
  const arma::uword n = X.n_rows;
  const arma::uword m = Time.n_elem;
  
  // precompute
  arma::vec exp_beta_x = ComputeExpBetaX(beta, X);
  arma::vec tailLam = tail_cumsum(lambda); // length m+1
  
  // --- First pass: for each subject i, compute the range of k covered by its interval (L_i, R_i] and the constant factor c_i,
  //             and accordingly accumulate colsum_w(k) = sum_i w(i,k) = sum_i lambda_k * c_i
  arma::vec colsum_w(m, arma::fill::zeros);
  
  // for the later second term of the score (sum_{i,k} w(i,k) X_i), keep (k_start, k_end) and c_i for each i
  std::vector<std::pair<int,int>> krange(n, {-1,-2}); // default: empty interval
  arma::vec ci(n, arma::fill::zeros);
  
  for(arma::uword i = 0; i < n; ++i){
    int id_i = static_cast<int>(simpdat(i,0)); // 1-based index
    double L_i = simpdat(i,1);
    double R_i = simpdat(i,2);
    
    if(std::isinf(R_i)){ // right-censored: w(i,k) = 0, not included in the E-step weights
      krange[i] = {-1,-2}; // empty
      ci(i) = 0.0;
      continue;
    }
    
    arma::uword idxL = first_leq_desc(Time, L_i);     // first t_k <= L_i
    arma::uword idxR = first_leq_desc(Time, R_i);     // first t_k <= R_i
    // the interval (L, R] corresponds to k in [idxR, idxL-1]
    int k_start = (idxR < m) ? static_cast<int>(idxR) : -1;
    int k_end   = (idxL > 0) ? static_cast<int>(idxL - 1) : -1;
    
    if(k_start < 0 || k_end < 0 || k_start > k_end){
      krange[i] = {-1,-2};
      ci(i) = 0.0;
      continue;
    }
    
    double e_bx = exp_beta_x(id_i - 1);
    double SL = e_bx * ( (idxL < m) ? tailLam(idxL) : 0.0 );   // ∑_{t<=L} λ * exp(bx)
    double SR = e_bx * ( (idxR < m) ? tailLam(idxR) : 0.0 );   // ∑_{t<=R} λ * exp(bx)
    double denom = 1.0 - std::exp(SL - SR);
    if(denom < 1e-12) denom = 1e-12; // numerical stability
    
    double c = e_bx / denom;
    ci(i) = c;
    krange[i] = {k_start, k_end};
    
    for(int k = k_start; k <= k_end; ++k){
      colsum_w((arma::uword)k) += lambda((arma::uword)k) * c; // contribution of w(i,k) to colsum
    }
  }
  
  // --- Second step: slide the risk set, accumulating InfMat and the first term of the score (related to the risk-set mean mu)
  arma::mat InfMat(p, p, arma::fill::zeros);
  arma::vec Score(p, arma::fill::zeros);
  
  double temp0 = 0.0;            // ∑_{j∈risk} exp(x_j^T β)
  arma::vec temp1(p, arma::fill::zeros);             // ∑ exp * x
  arma::mat temp2(p, p, arma::fill::zeros);          // ∑ exp * x x^T
  arma::uword idx_star = 0; // points to the newly added member of simpdat satisfying simpdat(j,3) >= t_k
  
  for(arma::uword k = 0; k < m; ++k){
    double t_k = Time(k);
    
    // add all subjects satisfying simpdat(idx_star,3) >= t_k to the risk set
    while(idx_star < n && simpdat(idx_star, 3) >= t_k){
      int id_j = static_cast<int>(simpdat(idx_star, 4)); // the 5th column (index 4) is the row number of X (1-based)
      arma::rowvec Xj = X.row(id_j - 1);
      double ebx = std::exp(arma::as_scalar(Xj * beta));
      temp0 += ebx;
      temp1 += ebx * arma::trans(Xj);
      temp2 += ebx * (arma::trans(Xj) * Xj);
      ++idx_star;
    }
    
    if(colsum_w(k) > 0 && temp0 > 0){
      arma::vec mu = temp1 / temp0; // weighted mean of x over the risk set
      arma::mat Sigma = (temp2 / temp0) - (mu * mu.t());
      
      InfMat += colsum_w(k) * Sigma;
      Score  += - colsum_w(k) * mu;  // first part of the score
    }
  }
  
  // --- Third step: second part of the score (sum_{i,k} w(i,k) X_i), accumulated directly without storing a per-k vector
  for(arma::uword i = 0; i < n; ++i){
    double c = ci(i);
    if(c <= 0.0) continue;
    int id = static_cast<int>(simpdat(i,0));
    arma::vec Xi = arma::trans(X.row(id - 1));
    auto pr = krange[i];
    for(int k = pr.first; k <= pr.second; ++k){
      Score += (lambda((arma::uword)k) * c) * Xi;
    }
  }
  
  // --- update beta
  arma::vec beta_new;
  // InfMat is theoretically symmetric positive definite, so use solve; add a small regularization if numerical issues occur
  arma::mat InfMat_reg = InfMat;
  // defensive regularization (can be reduced or removed as needed)
  InfMat_reg.diag() += 1e-8;
  beta_new = beta + arma::solve(InfMat_reg, Score, arma::solve_opts::likely_sympd);
  
  // --- update lambda (numerator uses colsum_w; denominator is the risk-set sum exp(x^T beta_new) under the new beta)
  arma::vec exp_beta_x_new = ComputeExpBetaX(beta_new, X);
  arma::vec denom_lambda(m, arma::fill::zeros);
  
  // slide the risk set again to obtain the denominator for each k
  temp0 = 0.0;
  idx_star = 0;
  for(arma::uword k = 0; k < m; ++k){
    double t_k = Time(k);
    while(idx_star < n && simpdat(idx_star, 3) >= t_k){
      int id_j = static_cast<int>(simpdat(idx_star, 4));
      temp0 += exp_beta_x_new(id_j - 1);
      ++idx_star;
    }
    denom_lambda(k) = std::max(temp0, 1e-16);
  }
  
  arma::vec lambda_new(m, arma::fill::zeros);
  for(arma::uword k = 0; k < m; ++k){
    lambda_new(k) = (denom_lambda(k) > 0.0) ? (colsum_w(k) / denom_lambda(k)) : 0.0;
  }
  
  MstepResult out;
  out.beta = beta_new;
  out.lambda = lambda_new;
  out.cumulative_lambda = CumulativeLambda(lambda_new);
  out.exp_beta_x = exp_beta_x_new;
  return out;
}

// ====== Main EM loop (reusing the memory-friendly M-step) ======
static EMResult EM(const arma::mat& X,
                   const arma::mat& simpdat,
                   double eps,
                   const arma::vec& Time){
  const arma::uword p = X.n_cols;
  const arma::uword m = Time.n_elem;
  
  arma::vec beta_hat(p, arma::fill::zeros);
  arma::vec lambda_hat(m, arma::fill::ones);
  lambda_hat /= (double)m;
  
  arma::vec cumulative_lambda_hat = CumulativeLambda(lambda_hat);
  arma::vec exp_beta_x_old = ComputeExpBetaX(beta_hat, X);
  
  double error = 100.0;
  int Iter = 0;
  const int MaxIter = 5000;
  const double delta = 0.01;
  
  while (error > eps && Iter < MaxIter) {
    arma::vec beta_old = beta_hat;
    arma::vec cumulative_lambda_old = cumulative_lambda_hat;
    
    MstepResult st = Mstep_memory_friendly(beta_hat, X, simpdat, Time, lambda_hat);
    
    beta_hat = st.beta;
    lambda_hat = st.lambda;
    cumulative_lambda_hat = st.cumulative_lambda;
    exp_beta_x_old = st.exp_beta_x;
    
    double beta_max = VectorAbsoluteDifferenceSafe(beta_hat, beta_old, delta);
    double cumlam_max = VectorAbsoluteDifferenceSafe(cumulative_lambda_hat, cumulative_lambda_old, delta);
    error = beta_max + cumlam_max;
    
    ++Iter;
    if(Iter == 1){
      Rcpp::Rcout << "\nAfter " << Iter  <<" iteration, beta=\n" << beta_hat << std::endl;
    }
    if(Iter % 50 == 1){
      Rcpp::Rcout << "Iter = " << Iter << ", Error = " << error << std::endl;
    }
  }
  
  if(error <= eps){
    Rcpp::Rcout << "Algorithm converged in " << Iter << " iterations, Error = " << error << std::endl;
  } else{
    Rcpp::Rcout << "Algorithm did not converge!" << std::endl;
  }
  
  EMResult res;
  res.beta = beta_hat;
  res.lambda = lambda_hat;
  res.cumulative_lambda = CumulativeLambda(lambda_hat);
  return res;
}


// [[Rcpp::export]]
Rcpp::List unireg_indep_EM_improved(arma::mat X, arma::mat simpdat, double eps, arma::vec Time){
  Rcpp::List ret;
  
  // EM
  EMResult res = EM(X, simpdat, eps, Time);
  
  
  ret["beta"] = res.beta;
  ret["lambda"] = res.lambda;
  ret["cumulative_lambda"] = res.cumulative_lambda;
  
  return ret;
}

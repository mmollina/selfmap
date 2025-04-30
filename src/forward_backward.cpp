// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

// Numerically stable log-sum-exp
static double log_sum_exp(const NumericVector & x) {
  double m = x[0];
  for (int i = 1; i < x.size(); ++i) if (x[i] > m) m = x[i];
  double s = 0.0;
  for (int i = 0; i < x.size(); ++i) s += std::exp(x[i] - m);
  return m + std::log(s);
}

// Initial state probabilities for F_n
static NumericVector initial_probs(int t) {
  double a = 0.5 - std::pow(2.0, -t);
  double b = std::pow(2.0, -(t - 1));
  return NumericVector::create(a, b, a);
}

// Compute the 3×3 log-transition matrix for one interval
static NumericMatrix compute_log_HFt(double r, int t) {
  double denom  = std::pow(2.0, t - 1) - 1.0;
  double half   = 0.5;
  double decay1 = std::pow(half - r, t);
  double decay2 = std::pow(half - r * (1.0 - r), t - 1);

  double H11 = (1.0/denom) * (
    std::pow(2.0, t - 1)/(1 + 2*r) - 1.0
  - std::pow(2.0, t - 1)*decay1/(1 + 2*r)
    + std::pow(2.0, t - 2)*decay2
  );
  double H12 = (1.0/denom) * (1.0 - std::pow(2.0, t - 1)*decay2);
  double H13 = (1.0/denom) * (
    std::pow(2.0, t)*r/(1 + 2*r) - 1.0
  + std::pow(2.0, t - 1)*decay1/(1 + 2*r)
    + std::pow(2.0, t - 2)*decay2
  );
  double H21 = half - std::pow(2.0, t - 2)*decay2;
  double H22 = std::pow(2.0, t - 1)*decay2;

  NumericMatrix M(3, 3);
  M(0,0) = H11;  M(0,1) = H12;  M(0,2) = H13;
  M(1,0) = H21;  M(1,1) = H22;  M(1,2) = H21;
  M(2,0) = H13;  M(2,1) = H12;  M(2,2) = H11;

  // avoid log(0)
  double eps = std::numeric_limits<double>::epsilon();
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      M(i,j) = std::log(std::max(M(i,j), eps));

  return M;
}

// [[Rcpp::export]]
List forward_backward_cpp(IntegerVector geno, int t, NumericVector r) {
  int T_len = geno.size();
  if (r.size() != T_len - 1) {
    stop("r must be a vector of length = length(geno) - 1");
  }

  // 1) Build list of log-transition matrices
  std::vector<NumericMatrix> logA;
  logA.reserve(T_len - 1);
  for (int k = 0; k < T_len - 1; ++k) {
    logA.push_back(compute_log_HFt(r[k], t));
  }

  // 2) Initial & emission probabilities
  NumericVector pi = initial_probs(t);
  NumericVector log_pi = wrap(log(pi));  // log of initial_probs

  NumericMatrix logB(T_len, 3);
  const double NEG_INF = -std::numeric_limits<double>::infinity();
  for (int tt = 0; tt < T_len; ++tt) {
    if (geno[tt] == NA_INTEGER) {
      logB(tt,0) = 0; logB(tt,1) = 0; logB(tt,2) = 0;
    } else if (geno[tt] == 0) {
      logB(tt,0) = 0; logB(tt,1) = NEG_INF; logB(tt,2) = NEG_INF;
    } else if (geno[tt] == 1) {
      logB(tt,0) = NEG_INF; logB(tt,1) = 0; logB(tt,2) = NEG_INF;
    } else {
      logB(tt,0) = NEG_INF; logB(tt,1) = NEG_INF; logB(tt,2) = 0;
    }
  }

  // 3) Forward pass
  NumericMatrix log_alpha(T_len, 3);
  for (int i = 0; i < T_len; ++i)
    for (int j = 0; j < 3; ++j)
      log_alpha(i,j) = NEG_INF;

  for (int j = 0; j < 3; ++j)
    log_alpha(0,j) = log_pi[j] + logB(0,j);

  for (int tt = 1; tt < T_len; ++tt) {
    NumericMatrix & A_prev = logA[tt-1];
    for (int j = 0; j < 3; ++j) {
      NumericVector temps(3);
      for (int i = 0; i < 3; ++i)
        temps[i] = log_alpha(tt-1,i) + A_prev(i,j);
      log_alpha(tt,j) = log_sum_exp(temps) + logB(tt,j);
    }
  }

  // log-likelihood
  NumericVector last_row = NumericVector::create(
    log_alpha(T_len-1, 0),
    log_alpha(T_len-1, 1),
    log_alpha(T_len-1, 2)
  );
  double loglik = log_sum_exp(last_row);

  // 4) Backward pass
  NumericMatrix log_beta(T_len, 3);
  for (int i = 0; i < T_len; ++i)
    for (int j = 0; j < 3; ++j)
      log_beta(i,j) = NEG_INF;

  for (int j = 0; j < 3; ++j)
    log_beta(T_len-1, j) = 0;

  for (int tt = T_len-2; tt >= 0; --tt) {
    NumericMatrix & A_curr = logA[tt];
    for (int i = 0; i < 3; ++i) {
      NumericVector temps(3);
      for (int j = 0; j < 3; ++j)
        temps[j] = A_curr(i,j) + logB(tt+1,j) + log_beta(tt+1,j);
      log_beta(tt,i) = log_sum_exp(temps);
    }
  }

  // 5) Posterior γ
  NumericMatrix log_gamma(T_len, 3), gamma(T_len, 3);
  for (int i = 0; i < T_len; ++i) {
    for (int j = 0; j < 3; ++j)
      log_gamma(i,j) = log_alpha(i,j) + log_beta(i,j);

    NumericVector row = NumericVector::create(
      log_gamma(i,0), log_gamma(i,1), log_gamma(i,2)
    );
    double denom = log_sum_exp(row);
    for (int j = 0; j < 3; ++j)
      gamma(i,j) = std::exp(log_gamma(i,j) - denom);
  }

  // 6) Joint-transition ξ
  int n_int = T_len - 1;
  NumericVector xi_vec(n_int * 3 * 3);
  xi_vec.attr("dim") = Dimension(n_int, 3, 3);

  for (int k = 0; k < n_int; ++k) {
    NumericMatrix & A_tt = logA[k];
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        int idx = k + n_int * (i + 3 * j);
        xi_vec[idx] = std::exp(
          log_alpha(k,i) +
            A_tt(i,j) +
            logB(k+1,j) +
            log_beta(k+1,j) -
            loglik
        );
      }
    }
  }

  return List::create(
    _["loglik"] = loglik,
    _["gamma"]  = gamma,
    _["xi"]     = xi_vec
  );
}

// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

List forward_backward_cpp(IntegerVector geno, int t, NumericVector r);

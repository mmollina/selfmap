// -----------------------------------------------------------------------------
// estimate_recombination.cpp   (revised, stateless ξ)
// -----------------------------------------------------------------------------
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(Rcpp)]]

#include <Rcpp.h>
#include "forward_backward.h"      // forward_backward_cpp()
// If forward_backward_cpp() is in the same translation unit you can omit this.
using namespace Rcpp;

// -----------------------------------------------------------------------------
// [[Rcpp::export]]
List estimate_recombination_cpp(DataFrame                 geno_df,
                                double                    tol       = 1e-5,
                                int                       max_iter  = 1000,
                                Nullable<NumericVector>   r_init    = R_NilValue)
{
  // -- basic checks -----------------------------------------------------------
  IntegerVector F_gen = geno_df["F_gen"];
  int n_ind = F_gen.size();

  for (int i = 0; i < n_ind; ++i) {
    if (F_gen[i] > 2) {
      warning("Some individuals have F_gen > 2. "
                "This experimental version assumes F2; results may be biased.");
      break;
    }
  }

  // marker names
  CharacterVector all_names = geno_df.names();
  std::vector<std::string> markers;
  for (auto & nm : all_names)
    if (std::string(nm) != "F_gen") markers.push_back(std::string(nm));

    int n_markers   = markers.size();
    if (n_markers < 2) stop("Need at least two marker columns");
    int n_intervals = n_markers - 1;

    // initial r
    NumericVector r(n_intervals);
    if (r_init.isNull())
      std::fill(r.begin(), r.end(), 0.1);
    else {
      NumericVector r0(r_init);
      if (r0.size() != n_intervals)
        stop("r_init must have length = number of intervals");
      r = r0;
    }

    NumericVector r_old(n_intervals, R_PosInf);
    NumericVector sum_odd(n_intervals);
    int iter = 0;

    // ————————————————————— EM LOOP ————————————————————————————————
    while (iter < max_iter) {

      // check convergence
      bool done = true;
      for (int k = 0; k < n_intervals; ++k)
        if (std::abs(r[k] - r_old[k]) > tol) { done = false; break; }
        if (done) break;

        r_old = clone(r);
        std::fill(sum_odd.begin(), sum_odd.end(), 0.0);

        // -------- E-step: compute expected # odd haplotypes for each interval ----
        for (int i = 0; i < n_ind; ++i) {

          // build genotype vector for individual i
          IntegerVector geno_i(n_markers);
          for (int k = 0; k < n_markers; ++k)
            geno_i[k] = IntegerVector(geno_df[ markers[k] ])[i];

          // run forward/backward with *current* r
          List fb   = forward_backward_cpp(geno_i, F_gen[i], r_old);
          NumericVector xi_arr = fb["xi"];
          IntegerVector dims   = xi_arr.attr("dim");   // [n_intervals,3,3]
          int di = dims[0];            // = n_intervals  (row stride)
          int dj = dims[1];            // = 3            (state stride)

          auto xi = [&](int k, int a, int b) {
            return xi_arr[k + di * (a + dj * b)];
          };

          for (int k = 0; k < n_intervals; ++k) {
            double p2 = 2 * r_old[k]*r_old[k] /
              (r_old[k]*r_old[k] + (1 - r_old[k])*(1 - r_old[k]));

            double agg = xi(k,0,1) + 2*xi(k,0,2) +
              xi(k,1,0) + p2*xi(k,1,1) +
              xi(k,1,2) + 2*xi(k,2,0) +
              xi(k,2,1);

            sum_odd[k] += agg;
          }
        }

        // divide by (2 gametes × n_ind) and update r
        for (int k = 0; k < n_intervals; ++k)
          r[k] = sum_odd[k] / (2.0 * n_ind);

        ++iter;
    }

    if (iter == max_iter)
      warning("Reached max_iter without full convergence");

    return List::create(_["r"] = r,
                        _["iterations"] = iter);
}

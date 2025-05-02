// -----------------------------------------------------------------------------
// estimate_recombination.cpp   (with delta-weighting)
// -----------------------------------------------------------------------------
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(Rcpp)]]

#include <Rcpp.h>
#include "forward_backward.h"
using namespace Rcpp;

// inline version of your R delta()
inline double delta_cpp(int t) {
  switch(t) {
  case 2: return 1.000;
  case 3: return 1.500;
  case 4: return 1.745;
  case 5: return 1.866;
  case 6: return 1.925;
  case 7: return 1.955;
  case 8: return 1.967;
  case 9: return 1.976;
  default: return 1.000;
  }
}

// [[Rcpp::export]]
List estimate_recombination_cpp(DataFrame               geno_df,
                                double                  tol      = 1e-5,
                                int                     max_iter = 1000,
                                Nullable<NumericVector> r_init   = R_NilValue)
{
  IntegerVector F_gen = geno_df["F_gen"];
  int n_ind = F_gen.size();

  // warn if not F2-ish
  for (int i = 0; i < n_ind; ++i) {
    if (F_gen[i] > 2) {
      warning("Some individuals have F_gen > 2; results may be biased.");
      break;
    }
  }

  // collect marker columns
  CharacterVector all_names = geno_df.names();
  std::vector<std::string> markers;
  for (auto & nm : all_names)
    if (std::string(nm) != "F_gen")
      markers.push_back(std::string(nm));

    int n_markers   = markers.size();
    if (n_markers < 2) stop("Need at least two marker columns");
    int n_intervals = n_markers - 1;

    // init r
    NumericVector r(n_intervals);
    if (r_init.isNull()) {
      std::fill(r.begin(), r.end(), 0.1);
    } else {
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
      // convergence?
      bool done = true;
      for (int k = 0; k < n_intervals; ++k) {
        if (std::abs(r[k] - r_old[k]) > tol) {
          done = false; break;
        }
      }
      if (done) break;

      r_old = clone(r);
      std::fill(sum_odd.begin(), sum_odd.end(), 0.0);

      // -------- E-step with delta weights -------------------------
      for (int i = 0; i < n_ind; ++i) {
        // genotype vector for individual i
        IntegerVector geno_i(n_markers);
        for (int k = 0; k < n_markers; ++k)
          geno_i[k] = IntegerVector(geno_df[ markers[k] ])[i];

        // run forward/backward
        List fb          = forward_backward_cpp(geno_i, F_gen[i], r_old);
        NumericVector xi = fb["xi"];
        IntegerVector dims = xi.attr("dim");   // [n_intervals,3,3]
        int di = dims[0], dj = dims[1];

        auto get_xi = [&](int k, int a, int b) {
          return xi[k + di * (a + dj * b)];
        };

        double d = delta_cpp(F_gen[i]);
        for (int k = 0; k < n_intervals; ++k) {
          // parity weight
          double p2 = 2 * r_old[k]*r_old[k] /
            (r_old[k]*r_old[k] + (1 - r_old[k])*(1 - r_old[k]));

          // aggregate xi slice
          double agg = get_xi(k,0,1) + 2*get_xi(k,0,2) +
            get_xi(k,1,0) + p2*get_xi(k,1,1) +
            get_xi(k,1,2) + 2*get_xi(k,2,0) +
            get_xi(k,2,1);

          // apply delta-scaling and F>9 correction
          double temp = agg / (2.0 * n_ind * d);
          if (F_gen[i] > 9)
            temp = temp / (2.0 * (1.0 - temp));

          sum_odd[k] += temp;
        }
      }

      // — M-step: r = sum_odd (already scaled inside the loop) —
      for (int k = 0; k < n_intervals; ++k) {
        r[k] = sum_odd[k];
      }

      ++iter;
    }

    if (iter == max_iter)
      warning("Reached max_iter without full convergence");

    return List::create(_["r"]          = r,
                        _["iterations"] = iter);
}

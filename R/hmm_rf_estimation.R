#' Estimate recombination fractions with an HMM–EM algorithm
#'
#' `estimate_recombination()` is a thin R-level wrapper around the
#' C++ engine `estimate_recombination_cpp()`.
#' The routine treats each chromosome as a three-state hidden-Markov model
#' (AA, AB, BB) whose transition matrix is a function of the unknown vector
#' **r** = (r\_1, …, r\_{m-1}).
#' An Expectation–Maximisation (EM) loop is run:
#'
#' * **E-step** – for every individual the C++ routine calls
#'   `forward_backward_cpp()` to obtain the posterior joint-transition
#'   probabilities \eqn{\xi\_{k,ij}} at the current **r**.
#' * **M-step** – the expected numbers of *odd* vs *even* recombinant
#'   gametes are updated and a closed-form parity estimator is used to
#'   refresh \eqn{r\_k}.
#'
#' The algorithm converges when the maximum change in **r** across
#' intervals is below `tol` or when `max_iter` iterations are reached.
#'
#' @section Important assumptions:
#' * The implementation is **optimised for F\eqn{_2} individuals**
#'   (one generation of selfing).  It will run on later generations
#'   (F\eqn{_3}, F\eqn{_6}, …) but the estimates can be biased.
#' * All marker columns must have exactly the same ordering for every
#'   individual and contain the discrete dosages 0 / 1 / 2 or `NA`.
#'
#' @param geno_df  A `data.frame` with one row per individual.
#'                 The first column must be named **`F_gen`** and give the
#'                 selfing generation number (e.g. `2` for F\eqn{_2}).
#'                 Every other column is interpreted as a biallelic marker
#'                 encoded 0, 1, 2 or `NA`.
#' @param tol      Convergence threshold for the maximum absolute change in
#'                 \eqn{r\_k} between successive EM iterations.
#'                 Default is `1e-5`.
#' @param max_iter Maximum number of EM iterations.  Default `1000`.
#' @param r_init   Optional numeric vector of length *(markers – 1)* providing
#'                 starting values for **r**.  If `NULL`, all intervals start
#'                 at 0.1.
#'
#' @return A list with two components
#' \describe{
#'   \item{`r`}{A numeric vector of length *(markers – 1)* with the
#'              estimated recombination fractions.}
#'   \item{`iterations`}{The number of EM iterations actually performed.}
#' }
#'
#' @examples
#' ## Toy example with simulated F2 data
#' set.seed(123)
#' sim <- SIMpoly::simulate_selfing_multi(n.mrk = 10,
#'                                        map.len = 50,
#'                                        n.ind  = 50,
#'                                        F.generations = 2)
#' est <- estimate_recombination(sim$geno$F_2)
#' est$r
#'
#' @seealso
#' * `forward_backward_cpp()` – C++ forward–backward used in the E-step
#' * `estimate_recombination_cpp()` – the underlying C++ implementation
#'
#' @export
estimate_recombination <- function(geno_df,
                                   tol      = 1e-5,
                                   max_iter = 1000,
                                   r_init   = NULL) {
  estimate_recombination_cpp(geno_df, tol, max_iter, r_init)
}

delta <- function(t,r){
  if(t == 2) return(1.0)
  else if(t == 3) return(1.502)
  else if(t == 4) return(1.746)
  else if(t == 5) return(1.873)
  else if(t == 6) return(1.936)
  else if(t > 6) return(1.0)
}

#' @export
estimate_recombination_R_version <- function(geno_df,
                                             tol      = 1e-5,
                                             max_iter = 1000,
                                             r_init   = NULL) {
  if (any(geno_df$F_gen > 2)) {
    warning(
      "Some individuals have F_gen > 2. ",
      "This experimental version of the algorithm assumes F2 (one selfing cycle), ",
      "so results for later generations (e.g., F3, F6) may be biased. ",
      "Use with caution."
    )
  }
  marker_cols <- setdiff(names(geno_df), "F_gen")
  n_markers   <- length(marker_cols)
  if (n_markers < 2) stop("Need at least two marker columns")
  n_intervals <- n_markers - 1

  if (is.null(r_init)) {
    r_init <- rep(0.1, n_intervals)
  }
  if (length(r_init) != n_intervals) {
    stop("r_init must have length = number of intervals")
  }

  n_ind   <- nrow(geno_df)
  # Precompute all xi's once per individual
  Xi_list <- lapply(seq_len(n_ind), function(i) {
    t_i    <- geno_df[["F_gen"]][i]
    geno_i <- as.numeric(geno_df[i, marker_cols])
    forward_backward(geno_i, t = t_i, r = r_init)$xi
  })

  r_old <- rep(Inf, n_intervals)
  r     <- r_init
  iter  <- 0

  while (any(abs(r - r_old) > tol) && iter < max_iter) {
    r_old <- r
    sum_odd <- numeric(n_intervals)

    for (k in seq_len(n_intervals)) {
      for (j in seq_len(n_ind)) {
        xi_k <- Xi_list[[j]][k, , ]
        d <- delta(geno_df[["F_gen"]][j])
        # parity weight for AA→AA etc.
        p2 <- 2 * r_old[k]^2 / (r_old[k]^2 + (1 - r_old[k])^2)
        temp <- ((
          xi_k[1, 2] + 2 * xi_k[1, 3] +
            xi_k[2, 1] + p2 * xi_k[2, 2] +
            xi_k[2, 3] + 2 * xi_k[3, 1] +
            xi_k[3, 2]
        ) / (2 * n_ind * d))
        if(geno_df[["F_gen"]][j] > 6)
          temp <-   temp / (2 * (1 - temp))
        sum_odd[k] <- sum_odd[k] + temp
      }
    }
    r <- sum_odd
    iter <- iter + 1
  }

  if (iter == max_iter) {
    warning("Reached max iterations without full convergence")
  }

  list(r = r, iterations = iter)
}

#' @title Recombination‐Fraction Estimation via EM
#' @description Given a data frame of genotypes and F_gen, iteratively estimate the vector of r's.
#' @param geno_df  Data.frame with columns “F_gen” plus marker columns (0/1/2/NA)
#' @param tol      Convergence tolerance (default 1e-5)
#' @param max_iter Max EM iterations (default 1e3)
#' @param r_init   Initial vector of length (markers − 1)
#' @return A list with `r` (estimated recombination fractions) and `iterations`
#' @export
estimate_recombination <- function(geno_df,
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
        # parity weight for AA→AA etc.
        p2 <- 2 * r_old[k]^2 / (r_old[k]^2 + (1 - r_old[k])^2)
        sum_odd[k] <- sum_odd[k] + (
          xi_k[1, 2] + 2 * xi_k[1, 3] +
            xi_k[2, 1] + p2 * xi_k[2, 2] +
            xi_k[2, 3] + 2 * xi_k[3, 1] +
            xi_k[3, 2]
        ) / (2 * n_ind)
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

#' @title Forward–Backward Algorithm
#' @description Compute log-likelihood, posterior state probabilities (γ), and joint-transition probabilities (ξ) for a three-state HMM over selfing generations.
#' @param geno Integer vector of observations (0, 1, 2 or NA).
#' @param t    Number of selfing generations (e.g., F2 → t=2).
#' @param r    Numeric vector of recombination fractions (length = length(geno) – 1).
#' @return A list with elements:
#'   \describe{
#'     \item{loglik}{Log-likelihood of the observed sequence.}
#'     \item{gamma}{T×3 matrix of posterior state probabilities at each marker.}
#'     \item{xi}{(T-1)×3×3 array of joint-transition probabilities between adjacent markers.}
#'   }
#' @export
#' @examples
#' # Example usage:
#' geno <- c(0, 1, NA, 2, 1, 0, 2)
#' t_val <- 2
#' r_vec <- rep(0.1, length(geno) - 1)
#' out <- forward_backward(geno, t = t_val, r = r_vec)
#'
#' # Inspect results:
#' out$loglik           # scalar log-likelihood
#' head(out$gamma)      # first few posterior probabilities
#' dim(out$xi)          # should be length(geno)-1, 3, 3
#' out$xi[5, , ]        # joint probabilities for interval 5→6
#'
#' # Plot posterior probabilities for each state over markers:
#' plot(out$gamma[, 1], type = "b", ylab = "P(state=1)", xlab = "Marker index", col = 2)
#' points(out$gamma[, 2], type = "b", col = 3)
#' points(out$gamma[, 3], type = "b", col = 4)

forward_backward <- function(geno, t, r) {
  T_len <- length(geno)
  if (length(r) != T_len - 1) {
    stop("r must be a vector of length = length(geno) - 1")
  }

  # Build array of log-transition matrices
  logA <- array(NA, dim = c(T_len - 1, 3, 3))
  for (k in seq_len(T_len - 1)) {
    logA[k, , ] <- compute_log_HFt(r[k], t)
  }

  # Initial and emission probabilities
  log_pi <- log(initial_probs(t))
  logB   <- matrix(0, nrow = T_len, ncol = 3)
  for (tt in seq_len(T_len)) {
    o <- geno[tt]
    if (is.na(o)) {
      logB[tt, ] <- 0
    } else if (o == 0) {
      logB[tt, ] <- diag(log_I0)
    } else if (o == 1) {
      logB[tt, ] <- diag(log_I1)
    } else {
      logB[tt, ] <- diag(log_I2)
    }
  }

  # Forward pass
  log_alpha <- matrix(-Inf, nrow = T_len, ncol = 3)
  log_alpha[1, ] <- log_pi + logB[1, ]
  for (tt in 2:T_len) {
    A_prev <- logA[tt - 1, , ]
    for (j in 1:3) {
      log_alpha[tt, j] <- log_sum_exp(log_alpha[tt - 1, ] + A_prev[, j]) +
        logB[tt, j]
    }
  }
  loglik <- log_sum_exp(log_alpha[T_len, ])

  # Backward pass
  log_beta <- matrix(-Inf, nrow = T_len, ncol = 3)
  log_beta[T_len, ] <- 0
  for (tt in (T_len - 1):1) {
    A_curr <- logA[tt, , ]
    for (i in 1:3) {
      log_beta[tt, i] <- log_sum_exp(
        A_curr[i, ] + logB[tt + 1, ] + log_beta[tt + 1, ]
      )
    }
  }

  # Posterior state probabilities (gamma)
  log_gamma <- log_alpha + log_beta
  gamma     <- exp(
    log_gamma - matrix(apply(log_gamma, 1, log_sum_exp),
                       nrow = T_len, ncol = 3, byrow = TRUE)
  )

  # Joint-transition probabilities (xi)
  xi <- array(0, dim = c(T_len - 1, 3, 3))
  for (tt in seq_len(T_len - 1)) {
    A_tt <- logA[tt, , ]
    for (i in 1:3) for (j in 1:3) {
      xi[tt, i, j] <- exp(
        log_alpha[tt, i] +
          A_tt[i, j] +
          logB[tt + 1, j] +
          log_beta[tt + 1, j] -
          loglik
      )
    }
  }

  list(loglik = loglik, gamma = gamma, xi = xi)
}

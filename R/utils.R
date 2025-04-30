#' @title Hidden‐Markov Model Utilities
#' @description Helper functions and constants for HMM‐based recombination estimation.
#' @keywords internal
#' @noRd

# Numerically stable log‐sum‐exp
log_sum_exp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

# Initial state probabilities for an Fₙ generation
initial_probs <- function(t) {
  a <- 1/2 - 1/2^t
  b <- 1/2^(t - 1)
  c(a, b, a)
}

# Log‐transition matrix (“logHFt”)
compute_log_HFt <- function(r, t) {
  denom  <- 2^(t - 1) - 1
  half   <- 1/2
  decay1 <- (half - r)^t
  decay2 <- (half - r * (1 - r))^(t - 1)

  H11 <- (1/denom) * (
    2^(t - 1)/(1 + 2*r) - 1 -
      2^(t - 1)*decay1/(1 + 2*r) +
      2^(t - 2)*decay2
  )
  H12 <- (1/denom) * (1 - 2^(t - 1)*decay2)
  H13 <- (1/denom) * (
    2^t * r/(1 + 2*r) - 1 +
      2^(t - 1)*decay1/(1 + 2*r) +
      2^(t - 2)*decay2
  )
  H21 <- half - 2^(t - 2)*decay2
  H22 <- 2^(t - 1)*decay2

  M <- matrix(c(
    H11, H12, H13,
    H21, H22, H21,
    H13, H12, H11
  ), nrow = 3, byrow = TRUE)

  log(pmax(M, .Machine$double.eps))
}

# Emission incidence matrices (in log‐space)
I0        <- diag(c(1, 0, 0))
I1        <- diag(c(0, 1, 0))
I2        <- diag(c(0, 0, 1))
I_missing <- I0 + I1 + I2

log_I0     <- {x <- log(I0);      x[!is.finite(x)] <- -Inf; x}
log_I1     <- {x <- log(I1);      x[!is.finite(x)] <- -Inf; x}
log_I2     <- {x <- log(I2);      x[!is.finite(x)] <- -Inf; x}
log_I_miss <- log(I_missing)  # all zeros

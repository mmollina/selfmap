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
  forward_backward_cpp(geno, t, r)
}

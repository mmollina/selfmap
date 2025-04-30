#' Whole-matrix HMM log-likelihood
#'
#' `log_likelihood_hmm()` evaluates the summed log-likelihood of a hidden-Markov
#' model for *all* individuals in a genotype matrix, given a fixed vector of
#' recombination fractions **r**.
#'
#' The function is most useful for **comparing candidate maps**: you supply the
#' same genotype matrix but different **r** vectors (e.g. Haldane distances,
#' EM-estimated distances, linkage-count distances) and see which one yields the
#' higher likelihood.
#'
#' @section Algorithm:
#' \enumerate{
#'   \item The first column (assumed to be \code{F_gen}) is removed.
#'   \item For each individual (row) the routine calls \code{forward_backward()}
#'         with the supplied selfing generation \code{t} and the common
#'         recombination-fraction vector \code{r}.
#'   \item The per-individual log-likelihoods are summed and returned.
#' }
#'
#' @param geno A matrix or data frame whose **first column** is \code{F_gen} and
#'   whose remaining columns are biallelic marker dosages encoded \code{0},
#'   \code{1}, \code{2}, or \code{NA}.
#' @param t    The selfing generation used for all individuals (e.g.\ \code{2}
#'   for F\subscript{2}).
#' @param r    **Recombination-fraction vector** of length
#'   \eqn{m-1}, where \eqn{m = ncol(geno) - 1} is the number of marker columns.
#'   Element \code{r[k]} is the recombination fraction between marker \code{k}
#'   and marker \code{k+1}.
#'
#' @return A single numeric value: the summed log-likelihood across
#'   individuals.
#'
#' @examples
#' ## Compare unbiased (F2) vs biased (F6) r-estimates
#' set.seed(2025)
#' res <- SIMpoly::simulate_selfing_multi(n.mrk = 20,
#'                                        map.len = 100,
#'                                        n.ind  = 300,
#'                                        F.generations = 6)
#'
#' ## “True” Haldane distances
#' r_true <- mappoly::mf_k(diff(res$map))
#'
#' ## Unbiased estimate from F2 only
#' r_F2 <- estimate_recombination(res$geno$F_2)$r
#'
#' l_true <- log_likelihood_hmm(res$geno$F_2, 2, r_true)
#' l_F2   <- log_likelihood_hmm(res$geno$F_2, 2, r_F2)
#' which.max(c(l_true, l_F2))  # 2 → EM estimate better than simulated distances
#'
#' ## Biased estimate when using only F6 individuals
#' r_F6 <- estimate_recombination(res$geno$F_6)$r
#' l_F6_true <- log_likelihood_hmm(res$geno$F_6, 6, r_true)
#' l_F6_est  <- log_likelihood_hmm(res$geno$F_6, 6, r_F6)
#' which.max(c(l_F6_true, l_F6_est))  # 1 → simulated distances better for F6
#'
#' @seealso `forward_backward()`, `estimate_recombination()`
#' @export
log_likelihood_hmm <- function(geno, t, r) {
  g <- geno[, -1, drop = FALSE]                      # remove F_gen column
  sum(apply(g, 1, function(x)
    forward_backward(x, t = t, r = r)$loglik))
}





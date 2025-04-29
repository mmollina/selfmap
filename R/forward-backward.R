# ─────────────────────────────────────────────────────────────────────────────
# 1. Helper functions & constants
# ─────────────────────────────────────────────────────────────────────────────

# your initial_probs
initial_probs <- function(t) {
  a <- 1/2 - 1/2^t
  b <- 1/2^(t - 1)
  c(a, b, a)
}

# your log-transition matrix
compute_log_HFt <- function(r, t) {
  denom <- 2^(t - 1) - 1
  half  <- 1/2
  decay1 <- (half - r)^t
  decay2 <- (half - r*(1 - r))^(t - 1)
  
  H11 <- (1/denom) * (
    (2^(t-1))/(1 + 2*r) - 1 -
      (2^(t-1)*decay1)/(1 + 2*r) +
      2^(t-2)*decay2
  )
  H12 <- (1/denom) * (1 - 2^(t - 1)*decay2)
  H13 <- (1/denom) * (
    (2^t * r)/(1 + 2*r) - 1 +
      (2^(t-1)*decay1)/(1 + 2*r) +
      2^(t-2)*decay2
  )
  H21 <- half - 2^(t-2)*decay2
  H22 <- 2^(t-1)*decay2
  
  M <- matrix(c(
    H11, H12, H13,
    H21, H22, H21,
    H13, H12, H11
  ), nrow = 3, byrow = TRUE)
  log(pmax(M, .Machine$double.eps))
}

# incidence (emission) matrices
I0 <- diag(c(1,0,0))
I1 <- diag(c(0,1,0))
I2 <- diag(c(0,0,1))
I_missing <- I0 + I1 + I2

# log‐incidences
log_I0 <- log(I0);      log_I0[!is.finite(log_I0)] <- -Inf
log_I1 <- log(I1);      log_I1[!is.finite(log_I1)] <- -Inf
log_I2 <- log(I2);      log_I2[!is.finite(log_I2)] <- -Inf
log_I_miss <- log(I_missing)  # this is all zeros

# numerically stable log‐sum‐exp
log_sum_exp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

# ─────────────────────────────────────────────────────────────────────────────
# Assumes you already have these on the search path:
#   initial_probs(t)
#   compute_log_HFt(r, t)
#   log_I0, log_I1, log_I2   # each is a 3×3 with log(1)/-Inf on the diag
#   log_sum_exp(v)           # numerically stable
# ─────────────────────────────────────────────────────────────────────────────
forward_backward <- function(geno, t, r) {
  T <- length(geno)
  if (length(r) != T-1) {
    stop("r must be a vector of length = length(geno) - 1")
  }
  
  # 1) build array of log-transition matrices:
  #    logA[k,,] is the 3×3 log-A for interval k→k+1
  logA <- array(NA, dim = c(T-1, 3, 3))
  for(k in 1:(T-1)) {
    logA[k,,] <- compute_log_HFt(r[k], t)
  }
  
  # 2) log-initial
  log_pi <- log(initial_probs(t))
  
  # 3) build T×3 log-emission matrix
  logB <- matrix(0, nrow=T, ncol=3)
  for(tt in seq_len(T)) {
    o <- geno[tt]
    if      (is.na(o)) logB[tt,] <- 0
    else if (o==0)     logB[tt,] <- c(log_I0[1,1], log_I0[2,2], log_I0[3,3])
    else if (o==1)     logB[tt,] <- c(log_I1[1,1], log_I1[2,2], log_I1[3,3])
    else               logB[tt,] <- c(log_I2[1,1], log_I2[2,2], log_I2[3,3])
  }
  
  # ───────────────────────────────────────────────────────────────────────────
  # Forward pass
  # ───────────────────────────────────────────────────────────────────────────
  log_alpha <- matrix(-Inf, nrow=T, ncol=3)
  log_alpha[1,] <- log_pi + logB[1,]
  
  for(tt in 2:T) {
    A_prev <- logA[tt-1,,]  # transitions from tt-1 → tt
    for(j in 1:3) {
      log_alpha[tt,j] <- log_sum_exp(log_alpha[tt-1,] + A_prev[,j]) + logB[tt,j]
    }
  }
  loglik <- log_sum_exp(log_alpha[T,])
  
  # ───────────────────────────────────────────────────────────────────────────
  # Backward pass
  # ───────────────────────────────────────────────────────────────────────────
  log_beta <- matrix(-Inf, nrow=T, ncol=3)
  log_beta[T,] <- 0
  
  for(tt in (T-1):1) {
    A_curr <- logA[tt,,]  # transitions from tt → tt+1
    for(i in 1:3) {
      log_beta[tt,i] <- log_sum_exp(
        A_curr[i,] + logB[tt+1,] + log_beta[tt+1,]
      )
    }
  }
  
  # ───────────────────────────────────────────────────────────────────────────
  # Posterior state probs γ
  # ───────────────────────────────────────────────────────────────────────────
  log_gamma <- log_alpha + log_beta
  gamma <- exp(
    log_gamma - matrix(apply(log_gamma,1,log_sum_exp), nrow=T, ncol=3)
  )
  
  # ───────────────────────────────────────────────────────────────────────────
  # Joint-transition probs ξ
  # ───────────────────────────────────────────────────────────────────────────
  xi <- array(0, dim = c(T-1,3,3))
  for(tt in 1:(T-1)) {
    A_tt <- logA[tt,,]
    for(i in 1:3) for(j in 1:3) {
      xi[tt,i,j] <- exp(
        log_alpha[tt,i] +
          A_tt[i,j] +
          logB[tt+1,j] +
          log_beta[tt+1,j] -
          loglik
      )
    }
  }
  
  list(
    loglik = loglik,
    gamma  = gamma,
    xi     = xi
  )
}


# ─────────────────────────────────────────────────────────────────────────────
# 4. Example usage
# ─────────────────────────────────────────────────────────────────────────────
# Suppose geno is one individual: e.g.
# geno <- c(0,1,NA,2,1,0,2)
# t    <- 2
# r    <- rep(0.1, length(geno)-1)
# out  <- forward_backward(geno, t, r)
# out$loglik          # scalar
# head(out$gamma)     # first few posterior state probs
# dim(out$xi)         # length(geno)-1, 3, 3
# out$xi[5,,]
# plot(out$gamma[,1], type="b", ylab="P(AB)", xlab="Marker index", col = 2)
# points(out$gamma[,2], type="b", ylab="P(AB)", xlab="Marker index", col = 3)
# points(out$gamma[,3], type="b", ylab="P(AB)", xlab="Marker index", col = 4)

estimate_recombination <- function(geno_df,
                                   tol      = 1e-5,
                                   max_iter = 1000,
                                   r_init   = NULL) {
  
  ## 2) Identify marker columns & set up initial r ------------------
  marker_cols <- setdiff(names(geno_df), "F_gen")
  n_markers   <- length(marker_cols)
  if (n_markers < 2) stop("Need at least two marker columns.")
  n_intervals <- n_markers - 1
  
  if (is.null(r_init)) {
    r_init <- rep(0.1, n_intervals)
  }
  if (length(r_init) != n_intervals) {
    stop("r_init must have length = number of intervals (markers - 1).")
  }
  
  n_ind <- nrow(geno_df)
  
  ## 3) Precompute all xi arrays with forward–backward ---------------
  Xi_list <- vector("list", n_ind)
  for (i in seq_len(n_ind)) {
    t_i    <- geno_df[["F_gen"]][i]
    geno_i <- as.numeric(geno_df[i, marker_cols])
    Xi_list[[i]] <- forward_backward(geno_i, t = t_i, r = r_init)$xi
  }
  
  ## 4) EM-style M-step (parity‐based) -----------------------------
  r_old <- rep(Inf, n_intervals)
  r     <- r_init
  iter  <- 0
  
  while (any(abs(r - r_old) > tol) && iter < max_iter) {
    r_old <- r
    
    for (k in seq_len(n_intervals)) {
      sum_odd <- 0
      
      for (j in seq_len(n_ind)) {
        xi_k   <- Xi_list[[j]][k, , ]   # 3×3 slice
        t_i    <- geno_df[["F_gen"]][j]
        #n_self <- t_i - 1               # selfing cycles (F₂ → 1, F₆ → 5, etc.)
        
        # parity probabilities for *one* gamete:
        #delta  <- (1 - 2 * r_old[k])^n_self
        #p_odd  <- (1 - delta) / 2
        #p_even <- 1 - p_odd
        
        # weight for AB→AB (0 or 2 odd haplotypes):
        #p2 <- p_odd^2 / (p_even^2 + p_odd^2)
        p2 <- 2 * r_old[k]^2 / (r_old[k]^2 + (1-r_old[k])^2)
        if(t_i == 2){
          # accumulate expected # odd haplotypes
          sum_odd <-
            sum_odd +
            (xi_k[1, 2] +        # AA→AB 
               2 * xi_k[1, 3] +    # AA→BB  
               xi_k[2, 1] +        # AB→AA  
               p2 * xi_k[2, 2] +   # AB→AB  
               xi_k[2, 3] +        # AB→BB  
               2 * xi_k[3, 1] +    # BB→AA  
               xi_k[3, 2]) / (2 * n_ind)          # BB→AB  
        } else if (t_i > 2) {
          sum_odd <-
            sum_odd +
            (xi_k[1, 2] +        # AA→AB 
               2 * xi_k[1, 3] +    # AA→BB  
               xi_k[2, 1] +        # AB→AA  
               p2 * xi_k[2, 2] +   # AB→AB  
               xi_k[2, 3] +        # AB→BB  
               2 * xi_k[3, 1] +    # BB→AA  
               xi_k[3, 2]) / (2 * n_ind)  # BB→AB 
          #sum_odd <- (sum_odd/(2*(1-sum_odd)))
        } else {
          stop("Should not get here!")
        }
      }
      # divide by 2 gametes × n_ind individuals to get recombination fraction
      r[k] <- sum_odd 
    }
    
    iter <- iter + 1
  }
  
  if (iter == max_iter) {
    warning("estimate_recombination: reached max_iter without full convergence")
  }
  
  list(r          = r,
       iterations = iter)
}


plot_simulation<-function(df){
  
  # Fit linear model
  fit <- lm(estimated ~ simulated, data = df)
  slope <- round(coef(fit)[2], 3)
  
  # Find common maximum for axes
  max_val <- max(df$simulated, df$estimated, na.rm = TRUE)
  breaks_seq <- seq(0, ceiling(max_val/10)*10, by = 10)
  
  df[df==0] <- NA
  
  # Plot
  p <- ggplot(df, aes(x = simulated, y = estimated)) +
    
    # Density-colored points
    geom_pointdensity(adjust = 1/2, size = 3) +
    scale_color_viridis_c(option = "plasma", direction = -1, name = "Density") +
    
    # 1:1 reference line
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", linewidth = 1) +
    
    # Regression line
    geom_smooth(method = "lm", se = FALSE, color = "#d62728", linewidth = 1.5) +
    
    # Annotate slope
    annotate("label",
             x = 0.7 * max_val,
             y = 0.1 * max_val,
             label = paste0("Slope = ", slope),
             size = 5,
             fill = "white",
             color = "#d62728",
             label.size = 0.4,
             fontface = "bold") +
    
    # Labels and scales
    labs(
      title = "Simulated vs Estimated Genetic Map Positions",
      subtitle = "Dashed line = perfect match (y = x)",
      x = "Simulated Position (cM)",
      y = "Estimated Position (cM)"
    ) +
    scale_x_continuous(breaks = breaks_seq, limits = c(0, max(breaks_seq))) +
    scale_y_continuous(breaks = breaks_seq, limits = c(0, max(breaks_seq))) +
    
    # Force square plot
    coord_fixed(ratio = 1) +
    
    # Fancier theme
    theme_minimal(base_size = 16) +
    theme(
      plot.title      = element_text(face = "bold", size = 18, hjust = 0.5),
      plot.subtitle   = element_text(size = 14, hjust = 0.5, color = "gray30"),
      axis.title.x    = element_text(face = "bold", size = 16),
      axis.title.y    = element_text(face = "bold", size = 16),
      axis.text       = element_text(size = 13),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
  
  # Print plot
  print(p)
  
  
}

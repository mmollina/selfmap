#!/usr/bin/env Rscript
#===============================================================================
# Title:   Map‐Estimation Simulations
# Author:  Marcelo Mollinari
# Date:    2025‐05‐03
# Purpose: Run three simulation scenarios, produce grids of
#          simulated vs. estimated map positions, and save results.
#===============================================================================

library(selfmap)
library(SIMpoly)
library(doParallel)
library(foreach)
library(ggplot2)

#-------------------------------------------------------------------------------
# Helper: run a single set of simulations for generations F2...F9
#-------------------------------------------------------------------------------
run_generations <- function(n.ind, n.rep, gens = 2:9, force_F2 = FALSE) {
  # Detect cores & register parallel backend
  cores <- parallel::detectCores() - 1
  cl    <- makeCluster(cores)
  registerDoParallel(cl)

  on.exit({
    stopCluster(cl)
    registerDoSEQ()
  }, add = TRUE)

  # Pre‐allocate list of data.frames
  out <- vector("list", length(gens))
  names(out) <- paste0("F_", gens)

  for (t in gens) {
    message("Running generation F", t)

    out[[paste0("F_", t)]] <- foreach(
      i        = seq_len(n.rep),
      .combine = rbind,
      .packages= c("SIMpoly", "selfmap")
    ) %dopar% {
      # simulate selfing
      sim <- simulate_selfing_multi(
        n.mrk         = 100,
        map.len       = 100,
        n.ind         = as.integer(n.ind),
        F.generations = t
      )
      geno <- sim$geno[[paste0("F_", t)]]

      # optionally force all as F2 for bias‐correction
      if (force_F2) geno$F_gen <- 2

      # estimate RF & build map
      rf      <- estimate_recombination(geno)$r
      est_map <- cumsum(c(0, imf_haldane_cM(rf)))

      data.frame(
        Marker    = colnames(geno)[-1],
        simulated = sim$map,
        estimated = round(est_map, 2),
        stringsAsFactors = FALSE
      )
    }
  }

  return(out)
}

#-------------------------------------------------------------------------------
# Scenario 1: Get Weights (force F2)
#-------------------------------------------------------------------------------
df_weights <- run_generations(
  n.ind    = 500,
  n.rep    = 5000,
  gens     = 2:9,
  force_F2 = TRUE
)

p1 <- plot_simulation_grid(df_weights, title_prefix = "", ncol = 4)
ggsave(
  "map_estimation_grid_weights.jpg",
  plot   = p1,
  width  = 3000,
  height = 2000,
  units  = "px",
  dpi    = 300
)
saveRDS(df_weights, "df_weights.rds")


#-------------------------------------------------------------------------------
# Scenario 2: Test Weights (no force)
#-------------------------------------------------------------------------------
df_test <- run_generations(
  n.ind    = 500,
  n.rep    = 5000,
  gens     = 2:9,
  force_F2 = FALSE
)

p2 <- plot_simulation_grid(df_test, title_prefix = "", ncol = 4)
ggsave(
  "map_estimation_grid_corrected_for_bias.jpg",
  plot   = p2,
  width  = 3000,
  height = 2000,
  units  = "px",
  dpi    = 300
)
saveRDS(df_test, "df_test_corrected_for_bias.rds")


#-------------------------------------------------------------------------------
# Scenario 3: Mixture families (F2 + F6 pooled)
#-------------------------------------------------------------------------------
run_mixture <- function(n.ind, n.rep) {
  cores <- parallel::detectCores() - 1
  cl    <- makeCluster(cores)
  registerDoParallel(cl)

  on.exit({
    stopCluster(cl)
    registerDoSEQ()
  }, add = TRUE)

  df_mix <- foreach(
    i        = seq_len(n.rep),
    .combine = rbind,
    .packages= c("SIMpoly", "selfmap")
  ) %dopar% {
    sim <- simulate_selfing_multi(
      n.mrk         = 100,
      map.len       = 100,
      n.ind         = as.integer(n.ind)
      # default F.generations = 2
    )
    geno <- rbind(sim$geno[["F_2"]], sim$geno[["F_6"]])

    rf      <- estimate_recombination(geno)$r
    est_map <- cumsum(c(0, imf_haldane_cM(rf)))

    data.frame(
      Marker    = colnames(geno)[-1],
      simulated = sim$map,
      estimated = round(est_map, 2),
      stringsAsFactors = FALSE
    )
  }

  return(df_mix)
}

df_mix <- run_mixture(n.ind = 500, n.rep = 10000)
p3     <- plot_simulation(df_mix)

ggsave(
  "map_estimation_grid_corrected_for_bias_F2-F6.jpg",
  plot   = p3,
  width  = 3000,
  height = 2000,
  units  = "px",
  dpi    = 300
)
saveRDS(df_mix, "df_mix_corrected_for_bias.rds")

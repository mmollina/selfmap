# selfmap

**selfmap** is an experimental R package for estimating recombination fractions in selfing populations using a hidden Markov model (HMM) and expectation-maximization (EM) algorithm. This implementation focuses on partially inbred populations (e.g., F₂) and is currently optimized for that use case.

> ⚠️ **Warning**: Although the code can handle individuals from generations beyond F₂ (e.g., F₆), the current algorithm assumes an F₂ structure. Results for later generations may be biased. Use with caution in such cases.

## Installation

To install `selfmap` and its dependencies, including the simulation tool [`SIMpoly`](https://github.com/mmollina/SIMpoly), use:

```r
# Install SIMpoly first if not already installed
devtools::install_github("mmollina/SIMpoly", dependencies = TRUE)

# Then install selfmap (once it's hosted)
devtools::install_github("mmollina/selfmap", dependencies = TRUE)
```

## Required Packages

```r
library(selfmap)
library(SIMpoly)
library(tidyverse)
```

---

## 🧪 Simulation 1: Estimating Map Using Only F₂ Individuals

This simulation estimates the genetic map from 200 independent F₂ populations.

```r
n.mrk <- 21
map.len <- 100
n.ind <- 100
n.rep <- 200
type <- "est_hmm"
df <- NULL

for (i in 1:n.rep) {
  cat("Sim: ", i, "\n")
  res <- simulate_selfing_multi(n.mrk = n.mrk, map.len = map.len, n.ind = n.ind, F.generations = 6, plot = FALSE)
  
  if (type == "est_hmm") {
    geno <- res$geno$F_2
    est_rf <- estimate_recombination(geno)
    est_map <- cumsum(c(0, imf_haldane_cM(est_rf$r)))
    
    df <- rbind(df, data.frame(
      Marker = colnames(res$geno$F_2)[-1],
      simulated = res$map,
      estimated = round(est_map, 2)
    ))
    
  } else if (type == "count_co") {
    est_rf_from_co <- round(res$count.CO$F_2 / (2 * n.ind), 3)
    est_pos_from_co <- cumsum(c(0, 100 * est_rf_from_co))
    
    df <- rbind(df, data.frame(
      Marker = colnames(res$geno$F_2)[-1],
      simulated = res$map,
      estimated = round(est_pos_from_co, 2)
    ))
  }
}

plot_simulation(df)
```

---

## 🧪 Simulation 2: Estimating Map Using Both F₂ and F₆ Individuals

This simulation tests performance using both early- and late-generation selfing (F₂ and F₆). It runs 100 replicates.

```r
n.mrk <- 20
map.len <- 100
n.ind <- 200
n.rep <- 100
type <- "est_hmm"
df <- NULL

for (i in 1:n.rep) {
  cat("Sim: ", i, "\n")
  res <- simulate_selfing_multi(n.mrk = n.mrk, map.len = map.len, n.ind = n.ind, F.generations = 6, plot = FALSE)
  
  if (type == "est_hmm") {
    geno <- rbind(res$geno$F_2, res$geno$F_6)
    est_rf <- estimate_recombination(geno)
    
    # Optional: adjust for RIL mapping (not currently active)
    # adj.ril <- function(r.obs) r.obs / (2 * (1 - r.obs))
    # rx <- adj.ril(est_rf$r)
    
    est_map <- cumsum(c(0, imf_haldane_cM(est_rf$r)))
    
    df <- rbind(df, data.frame(
      Marker = colnames(res$geno$F_6)[-1],
      simulated = res$map,
      estimated = round(est_map, 2)
    ))
    
  } else if (type == "count_co") {
    est_rf_from_co <- round(res$count.CO$F_6 / (2 * n.ind), 3)
    est_pos_from_co <- cumsum(c(0, 100 * est_rf_from_co))
    
    df <- rbind(df, data.frame(
      Marker = colnames(res$geno$F_6)[-1],
      simulated = res$map,
      estimated = round(est_pos_from_co, 2)
    ))
  }
}

plot_simulation(df)
```

---

## 📈 Plotting Results

The function `plot_simulation(df)` creates a scatter plot of simulated vs. estimated genetic positions, colored by point density and annotated with the slope of the fitted regression line. This allows you to visually assess the accuracy and bias in your recombination estimates.

---

## License

MIT © Marcelo Mollinari  
For questions or contributions, please [submit an issue](https://github.com/mmollina/selfmap/issues) or a pull request.
```

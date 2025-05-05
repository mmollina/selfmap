#' Plot simulated vs. estimated map positions
#'
#' @param df A data frame with numeric columns `simulated` and `estimated`
#'
#' @importFrom stats lm coef
#' @importFrom ggplot2 ggplot aes geom_abline geom_smooth annotate labs
#'   scale_x_continuous scale_y_continuous coord_fixed theme_minimal theme
#'   scale_color_viridis_c
#' @importFrom ggpointdensity geom_pointdensity
#'
#' @export
plot_simulation <- function(df) {
  fit     <- lm(estimated ~ simulated, data = df)
  slope   <- round(coef(fit)[2], 3)
  max_val <- max(c(df$simulated, df$estimated), na.rm = TRUE)
  breaks  <- seq(0, ceiling(max_val / 10) * 10, by = 10)

  df[df == 0] <- NA

  ggplot(df, aes(simulated, estimated)) +
    geom_pointdensity(adjust = 0.5, size = 3) +
    scale_color_viridis_c(option = "plasma", direction = -1, name = "Density") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    geom_smooth(method = "lm", se = FALSE) +
    annotate(
      "label",
      x     = 0.7 * max_val,
      y     = 0.1 * max_val,
      label = paste0("Slope = ", slope),
      fontface = "bold"
    ) +
    labs(
      title    = "Simulated vs Estimated Genetic Map Positions",
      subtitle = "Dashed line = y = x",
      x        = "Simulated (cM)",
      y        = "Estimated (cM)"
    ) +
    scale_x_continuous(breaks = breaks, limits = c(0, max(breaks))) +
    scale_y_continuous(breaks = breaks, limits = c(0, max(breaks))) +
    coord_fixed() +
    theme_minimal(base_size = 16) +
    theme(legend.position = "right")
}

#' Plot a grid of simulated vs. estimated map positions
#'
#' @param df_list A named list of data frames. Each element must have numeric
#'   columns `simulated` and `estimated`.
#' @param title_prefix Optional string to prepend to each panel title.
#' @param ncol Number of columns in the output grid (default: 4).
#' @return A patchwork object combining all panels.
#'
#' @importFrom stats lm coef
#' @importFrom ggplot2 ggplot aes geom_abline geom_smooth annotate labs
#'   scale_x_continuous scale_y_continuous coord_fixed theme_minimal theme
#'   element_text
#' @importFrom ggpointdensity geom_pointdensity
#' @importFrom patchwork wrap_plots
#'
#' @export
plot_simulation_grid <- function(df_list, title_prefix = "", ncol = 4) {
  plot_one <- function(data, name) {
    fit     <- lm(estimated ~ simulated, data = data)
    slope   <- round(coef(fit)[2], 3)
    max_val <- max(c(data$simulated, data$estimated), na.rm = TRUE)
    breaks  <- seq(0, ceiling(max_val / 10) * 10, by = 10)

    data[data == 0] <- NA

    ggplot(data, aes(simulated, estimated)) +
      geom_pointdensity(adjust = 0.5, size = 2) +
      scale_color_viridis_c(option = "plasma", direction = -1, name = "Density") +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      geom_smooth(method = "lm", se = FALSE, color = "blue", size = 0.8) +
      annotate(
        "label",
        x     = 0.7 * max_val,
        y     = 0.1 * max_val,
        label = paste0("Slope = ", slope),
        fontface = "bold",
        size    = 2
      ) +
      labs(
        title = paste0(title_prefix, name),
        x     = "Simulated (cM)",
        y     = "Estimated (cM)"
      ) +
      scale_x_continuous(breaks = breaks, limits = c(0, max(breaks))) +
      scale_y_continuous(breaks = breaks, limits = c(0, max(breaks))) +
      coord_fixed() +
      theme_minimal(base_size = 10) +
      theme(
        legend.position = "none",
        axis.text.x    = element_text(angle = 90, hjust = 1, vjust = 0.5)
      )
  }

  plots <- mapply(
    plot_one,
    df_list,
    names(df_list),
    SIMPLIFY = FALSE,
    USE.NAMES = FALSE
  )

  wrap_plots(plots, ncol = ncol)
}
#' Compare two sets of simulated vs. estimated maps in a grid
#'
#' @param df_test_list    Named list of data frames (each with `simulated` and `estimated`) for your test set.
#' @param df_weights_list Named list of data frames (same names) for your weighted set.
#' @param title_prefix    Optional string to prepend to each panel title.
#' @param ncol            Number of columns in the output grid (default: 4).
#' @return A patchwork object combining all comparison panels.
#'
#' @importFrom stats lm coef
#' @importFrom ggplot2 ggplot aes geom_abline geom_smooth annotate labs
#'   scale_x_continuous scale_y_continuous coord_fixed theme_minimal theme
#'   element_text
#' @importFrom ggpointdensity geom_pointdensity
#' @importFrom patchwork wrap_plots
#' @importFrom viridis plasma magma
#' @export
plot_simulation_compare_grid <- function(df_test_list,
                                         df_weights_list,
                                         title_prefix = "",
                                         ncol = 4) {
  stopifnot(
    identical(names(df_test_list), names(df_weights_list)),
    length(df_test_list) > 0
  )

  plot_one <- function(test_df, wt_df, name) {
    # fit both regressions
    fit_test <- lm(estimated ~ simulated, data = test_df)
    fit_wt   <- lm(estimated ~ simulated, data = wt_df)
    s1 <- round(coef(fit_test)[2], 3)
    s2 <- round(coef(fit_wt)[2],   3)

    # determine axis limits & breaks
    all_vals <- c(test_df$simulated, test_df$estimated,
                  wt_df$simulated,   wt_df$estimated)
    maxv     <- max(all_vals, na.rm = TRUE)
    brks     <- seq(0, ceiling(maxv/10)*10, by = 10)

    # mask zeros
    test_df[test_df == 0] <- NA
    wt_df[wt_df == 0]     <- NA

    ggplot() +
      # test points in plasma palette
      geom_pointdensity(
        data     = test_df,
        mapping  = aes(simulated, estimated),
        gradient = viridis::plasma(100),
        adjust   = 0.5,
        size     = 2
      ) +
      # weights points in magma palette
      geom_pointdensity(
        data     = wt_df,
        mapping  = aes(simulated, estimated),
        gradient = viridis::magma(100),
        adjust   = 0.5,
        size     = 2
      ) +
      # reference diagonal
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      # regression lines
      geom_smooth(
        data   = test_df,
        aes(simulated, estimated),
        method = "lm",
        se     = FALSE,
        color  = "#1b9e77", # a greenish
        size   = 0.8
      ) +
      geom_smooth(
        data   = wt_df,
        aes(simulated, estimated),
        method = "lm",
        se     = FALSE,
        color  = "#d95f02", # an orange
        size   = 0.8
      ) +
      # annotate both slopes
      annotate(
        "label",
        x     = 0.65 * maxv,
        y     = 0.15 * maxv,
        label = paste0("test slope = ", s1),
        color = "#1b9e77",
        size  = 2
      ) +
      annotate(
        "label",
        x     = 0.65 * maxv,
        y     = 0.08 * maxv,
        label = paste0("wt slope   = ", s2),
        color = "#d95f02",
        size  = 2
      ) +
      labs(
        title = paste0(title_prefix, name),
        x     = "Simulated (cM)",
        y     = "Estimated (cM)"
      ) +
      scale_x_continuous(breaks = brks, limits = c(0, max(brks))) +
      scale_y_continuous(breaks = brks, limits = c(0, max(brks))) +
      coord_fixed() +
      theme_minimal(base_size = 10) +
      theme(
        legend.position = "none",
        axis.text.x     = element_text(angle = 90, hjust = 1, vjust = 0.5)
      )
  }

  # build one plot per element
  plots <- mapply(
    plot_one,
    df_test_list,
    df_weights_list,
    names(df_test_list),
    SIMPLIFY = FALSE,
    USE.NAMES = FALSE
  )

  # combine with patchwork
  patchwork::wrap_plots(plots, ncol = ncol)
}

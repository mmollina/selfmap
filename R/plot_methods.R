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
  max_val <- max(df$simulated, df$estimated, na.rm = TRUE)
  breaks  <- seq(0, ceiling(max_val / 10) * 10, by = 10)

  df[df == 0] <- NA

  ggplot(df, aes(simulated, estimated)) +
    geom_pointdensity(adjust = 1/2, size = 3) +
    scale_color_viridis_c(option = "plasma", direction = -1, name = "Density") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    geom_smooth(method = "lm", se = FALSE) +
    annotate("label", x = 0.7 * max_val, y = 0.1 * max_val,
             label = paste0("Slope = ", slope),
             fontface = "bold") +
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

# ==============================================================================
# Diversity.R
# Hill number diversity estimation and extrapolation curves
# Study: Dung beetle diversity across burning regimes
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load required libraries
# ------------------------------------------------------------------------------
library(tidyverse)  # Data manipulation and visualization
library(iNEXT)      # Interpolation and extrapolation of Hill numbers
library(ggpubr)     # Publication-ready figure assembly
library(cowplot)    # Plot composition utilities

# ------------------------------------------------------------------------------
# 2. Define theme and color palette
# ------------------------------------------------------------------------------

# Shared plot theme (serif font, clean axes, subtle grid)
base_theme <- theme_classic(base_family = "serif", base_size = 18) +
  theme(
    axis.text        = element_text(color = "black"),
    axis.title       = element_text(color = "black", size = 20, face = "bold"),
    axis.line        = element_line(color = "black"),
    axis.title.y     = element_text(margin = margin(r = 10)),
    axis.title.x     = element_text(margin = margin(t = 10)),
    panel.grid.major = element_line(color = "gray95", linetype = "dashed")
  )

# Color palette: Null = light, Low = medium, High = intense
fire_colors <- c("orangered", "orange", "lightgoldenrod1")

# ------------------------------------------------------------------------------
# 3. Load and clean data
# ------------------------------------------------------------------------------
dungbeetles <- read.csv("Diversity/Diversity.csv", sep = ";")
dungbeetles <- dungbeetles[, -1]  # Remove row index column

# ------------------------------------------------------------------------------
# 4. Compute diversity estimates (iNEXT)
# ------------------------------------------------------------------------------
q_values    <- c(0, 1, 2)   # Diversity orders: richness, Shannon, Simpson
div_outputs <- map(q_values, ~ iNEXT(dungbeetles, q = .x,
                                     datatype = "abundance", endpoint = 600))

# Asymptotic diversity estimators
estimates_raw <- iNEXT(dungbeetles, datatype = "abundance")$AsyEst |>
  as.data.frame()

# Prepare estimator table for plotting
div_estimates <- estimates_raw |>
  rename(
    OriginalAssemblage = Assemblage,
    SE                 = `s.e.`      # Rename standard error for clarity
  ) |>
  rownames_to_column("Assemblage") |>
  mutate(
    Order        = rep(paste0("q", q_values), times = nrow(estimates_raw) / 3),
    OrderNumeric = rep(1:3,                   times = nrow(estimates_raw) / 3)
  )

# ------------------------------------------------------------------------------
# 5. Plot: diversity estimators across orders
# ------------------------------------------------------------------------------
est_plot <- ggplot(div_estimates,
                   aes(x = OrderNumeric, y = Estimator,
                       color = OriginalAssemblage)) +
  geom_point() +
  geom_line(linewidth = 0.7) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL),
                linetype = "dotted", linewidth = 0.6) +
  scale_x_continuous(breaks = 1:3, labels = paste0("q", q_values)) +
  scale_colour_manual(values = fire_colors) +
  labs(x = "Diversity orders", y = "Estimator", color = "Burning regimes") +
  base_theme

# ------------------------------------------------------------------------------
# 6. Plot: rarefaction/extrapolation curves per diversity order
# ------------------------------------------------------------------------------
curve_plots <- map2(div_outputs, q_values, function(output, q) {

  p <- ggiNEXT(output, type = 1) +
    scale_colour_manual(name = "Burning regime", values = fire_colors) +
    scale_fill_manual(name   = "Burning regime", values = fire_colors) +
    base_theme +
    theme(axis.title.y = element_blank()) +
    labs(fill = "Burning regime", color = "Burning regime")

  y_lim <- ggplot_build(p)$layout$panel_params[[1]]$y.range[2]
  label <- paste0(c("\u{00B0}", "\u{00B9}", "\u{00B2}")[q + 1], "D")

  # Remove x-axis labels for upper panels
  if (q < 2) {
    p <- p + theme(axis.text.x  = element_blank(),
                   axis.title.x = element_blank())
  }

  p + annotate("text", x = 620, y = y_lim * 0.97,
               label = label, size = 7, family = "serif")
})

# ------------------------------------------------------------------------------
# 7. Assemble and export final figure
# ------------------------------------------------------------------------------
combined_curves <- ggarrange(
  plotlist     = curve_plots,
  ncol         = 1,
  nrow         = 3,
  align        = "v",
  common.legend = TRUE,
  legend       = "top"
)

final_plot <- annotate_figure(
  combined_curves,
  left = text_grob("Species diversity", rot = 90, vjust = 0.7,
                   size = 20, face = "bold", family = "serif")
)

ggsave("figures/Diversity.tiff", final_plot,
       width = 9, height = 11, dpi = 600, bg = "white")

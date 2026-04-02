# ==============================================================================
# RankAbundance.R
# Rank-abundance distribution (RAD) curves
# Study: Dung beetle species rank-abundance across burning regimes
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load helper functions and libraries
# ------------------------------------------------------------------------------
source("RangeAbundance/SpecDist.R")
source("RangeAbundance/DetAbu.R")
source("RangeAbundance/UndAbu.R")

library(iNEXT)    # Diversity estimation (loaded for consistency)
library(Rcpp)     # C++ integration
library(ggplot2)  # Visualization
library(ggrepel)  # Non-overlapping text labels
library(stringr)  # String manipulation utilities
library(dplyr)    # Data manipulation
library(readxl)   # Excel reading (if needed)
library(ggpubr)   # Figure annotation utilities

# ------------------------------------------------------------------------------
# 2. Load and aggregate data by burning regime
# ------------------------------------------------------------------------------
ra <- read.csv("RangeAbundance/RangeAbundance.csv", sep = ";")

# Sum abundances across sites within each burning regime
ra2 <- ra |>
  group_by(Species) |>
  summarise(across(c(Null, Low, High), sum))

# ------------------------------------------------------------------------------
# 3. Compute rank-abundance distributions (detected species only)
# ------------------------------------------------------------------------------

# Null regime
out1 <- SpecDist(ra2$Null, "abundance") %>% subset(probability > 0)
out1$number      <- seq_len(nrow(out1))
out1$sp          <- as.numeric(rownames(out1))
out1$Species     <- ra2$Species[out1$sp]
out1$BurningRegimen <- "Null"

# Low regime
out2 <- SpecDist(ra2$Low, "abundance") %>% subset(probability > 0)
out2$number      <- seq_len(nrow(out2))
out2$sp          <- as.numeric(rownames(out2))
out2$Species     <- ra2$Species[out2$sp]
out2$BurningRegimen <- "Low"

# High regime
out3 <- SpecDist(ra2$High, "abundance") %>% subset(probability > 0)
out3$number      <- seq_len(nrow(out3))
out3$sp          <- as.numeric(rownames(out3))
out3$Species     <- ra2$Species[out3$sp]
out3$BurningRegimen <- "High"

# ------------------------------------------------------------------------------
# 4. Combine and filter to detected species only
# ------------------------------------------------------------------------------
outtotal <- rbind(out1, out2, out3) %>%
  subset(method == "detected") %>%
  mutate(BurningRegimen = factor(BurningRegimen, levels = c("Null", "Low", "High")))

# ------------------------------------------------------------------------------
# 5. Define custom theme
# ------------------------------------------------------------------------------
theme_custom <- theme_minimal(base_family = "serif") +
  theme(
    axis.title       = element_text(size = 20, face = "bold"),
    strip.text.x     = element_blank(),
    axis.text        = element_text(size = 18),
    legend.text      = element_text(size = 18),
    legend.title     = element_text(size = 20, face = "bold"),
    legend.position  = "top",
    axis.line        = element_line(colour = "black"),
    panel.grid       = element_blank(),
    axis.title.x     = element_text(margin = margin(t = 20)),
    axis.title.y     = element_blank(),
    panel.grid.major = element_line(color = "gray95", linetype = "dashed")
  )

# ------------------------------------------------------------------------------
# 6. Plot rank-abundance curves
# ------------------------------------------------------------------------------
p <- ggplot(outtotal, aes(x = number, y = probability,
                          color = BurningRegimen, shape = BurningRegimen)) +
  geom_point(size = 4) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("lightgoldenrod2", "orange", "orangered")) +

  # Species name labels (italic, angled)
  geom_text(aes(label = Species), size = 5, family = "serif", color = "black",
            angle = 45, hjust = -0.15, vjust = -0.1, fontface = "bold.italic") +

  scale_x_continuous("Species rank") +
  ylim(0, 0.6) +
  xlim(0, 19) +
  facet_grid(cols = vars(BurningRegimen)) +
  labs(col = "Burning regime", shape = "Burning regime") +
  xlab("Species rank") +
  theme_custom

# Add shared y-axis label
p <- annotate_figure(p,
                     left = text_grob("Proportional abundance", rot = 90,
                                      vjust = 0.3, size = 20, face = "bold",
                                      family = "serif"))

# ------------------------------------------------------------------------------
# 7. Export figure
# ------------------------------------------------------------------------------
ggsave("figures/RangeAbundance.tiff",
       width = 15, height = 9, dpi = 600, bg = "white")

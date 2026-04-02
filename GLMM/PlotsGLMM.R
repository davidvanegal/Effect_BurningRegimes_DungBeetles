# ==============================================================================
# PlotsGLMM.R
# Publication figures — GLMM results
# Study: Dung beetle abundance, diversity, removal, and functional groups
#        across burning regimes
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load required libraries
# ------------------------------------------------------------------------------

# Modeling & statistics
library(lme4)           # Mixed-effects models
library(nlme)           # Linear and nonlinear mixed-effects models
library(multcomp)       # Multiple comparisons (Tukey, etc.)
library(emmeans)        # Estimated marginal means
library(ggeffects)      # Marginal effects visualization
library(lmtest)         # Regression model testing
library(MASS)           # Statistical functions (glm.nb, etc.)
library(car)            # Companion to Applied Regression

# Data manipulation
library(tidyverse)      # ggplot2, dplyr, tidyr, readr, purrr, ...
library(data.table)     # High-performance data handling
library(psych)          # Descriptive statistics utilities
library(lazyWeave)      # LaTeX wrappers for tables
library(PupillometryR)  # Raincloud plot utilities

# Visualization
library(ggpubr)         # Publication-ready figure assembly
library(cowplot)        # Plot composition
library(ggridges)       # Ridgeline plots
library(gg.layers)      # Additional ggplot2 layers (geom_boxplot2)
library(vcd)            # Categorical data visualization
library(lattice)        # Trellis graphics

# ------------------------------------------------------------------------------
# 2. Load and prepare data
# ------------------------------------------------------------------------------
GLMMDB <- read.csv("GLMM/Data.csv", sep = ";")
GLMMDB$BurningRegime <- factor(GLMMDB$BurningRegime,
                               levels = c("Null", "Low", "High"))

# ------------------------------------------------------------------------------
# 3. Define shared plot themes
# ------------------------------------------------------------------------------

# Primary theme: classic with serif font
theme_serif <- theme_classic(base_family = "serif") +
  theme(
    axis.text        = element_text(size = 18),
    axis.title       = element_text(size = 20, face = "bold"),
    axis.title.y     = element_text(margin = margin(r = 10)),
    axis.title.x     = element_text(margin = margin(t = 10)),
    panel.grid.major = element_line(color = "gray95", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(color = "black"),
    panel.background = element_blank(),
    legend.position  = "none",
    legend.text      = element_text(size = 18),
    legend.title     = element_text(size = 20, face = "bold")
  )

# Secondary theme: no y-axis title (for panels in composite figures)
theme_serif2 <- theme_serif +
  theme(axis.title.y = element_blank())

# Shared fill palette
fill_vals <- c("lightgoldenrod1", "orange", "orangered")

# ==============================================================================
# FIGURE 1: Abundance and Hill-number diversity
# ==============================================================================

# ------------------------------------------------------------------------------
# Panel A: Total abundance by burning regime
# ------------------------------------------------------------------------------
plotDB <- ggplot(GLMMDB, aes(x = BurningRegime, y = Abundance,
                             fill = BurningRegime)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.6) +
  scale_fill_manual(values = fill_vals) +
  scale_y_continuous("Number of individuals",
                     breaks = seq(0, 200, 40), limits = c(0, 200)) +
  # Post-hoc letter annotations
  annotate("text", label = "a", size = 6, x = 1, y = 110, family = "serif") +
  annotate("text", label = "a", size = 6, x = 2, y =  60, family = "serif") +
  annotate("text", label = "a", size = 6, x = 3, y = 140, family = "serif") +
  xlab("Burning regime") +
  theme_serif

# ------------------------------------------------------------------------------
# Panel B: Hill-number diversity (q = 0, 1, 2)
# ------------------------------------------------------------------------------
qplot <- GLMMDB |>
  pivot_longer(cols = c(q0, q1, q2), names_to = "q", values_to = "Total") |>
  mutate(q = factor(q, levels = c("q0", "q1", "q2"),
                    labels = c("0", "1", "2"))) |>
  ggplot(aes(x = q, y = Total, fill = BurningRegime)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.6) +
  scale_fill_manual(values = fill_vals) +
  scale_y_continuous("Species diversity",
                     breaks = seq(0, 10, 2), limits = c(0, 11)) +
  labs(x = "Diversity order (q)", fill = "Burning regime") +
  # Post-hoc letter annotations (grouped by q and regime)
  annotate("text", label = "a", size = 6, x = 0.8, y = 11.0, family = "serif") +
  annotate("text", label = "a", size = 6, x = 1.0, y =  7.0, family = "serif") +
  annotate("text", label = "a", size = 6, x = 1.2, y = 10.0, family = "serif") +
  annotate("text", label = "a", size = 6, x = 1.8, y =  6.8, family = "serif") +
  annotate("text", label = "a", size = 6, x = 2.0, y =  6.8, family = "serif") +
  annotate("text", label = "a", size = 6, x = 2.2, y =  6.8, family = "serif") +
  annotate("text", label = "a", size = 6, x = 2.8, y =  4.0, family = "serif") +
  annotate("text", label = "a", size = 6, x = 3.0, y =  6.0, family = "serif") +
  annotate("text", label = "a", size = 6, x = 3.2, y =  5.8, family = "serif") +
  theme_serif +
  theme(legend.position = "top")

# Combine panels and export
figura_plot1 <- ggarrange(plotDB, qplot, ncol = 2, nrow = 1, align = "v")
ggsave("figures/AbundDivers.svg", figura_plot1,
       width = 16, height = 9, dpi = 600, bg = "white")

# ==============================================================================
# FIGURE 2: Dung removal
# ==============================================================================

# ------------------------------------------------------------------------------
# Panel A: Community-level dung removal
# ------------------------------------------------------------------------------
plotReT <- ggplot(GLMMDB, aes(x = BurningRegime, y = RemovalCommunity,
                              fill = BurningRegime)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.6) +
  scale_y_continuous(breaks = seq(0, 350, by = 50), limits = c(0, 350)) +
  scale_fill_manual(values = fill_vals) +
  # Post-hoc letter annotations
  annotate("text", label = "a", size = 6, x = 1, y = 130, family = "serif") +
  annotate("text", label = "a", size = 6, x = 2, y = 210, family = "serif") +
  annotate("text", label = "a", size = 6, x = 3, y = 350, family = "serif") +
  labs(x = "Burning regime", y = "Dung removal by dung beetles (g)") +
  theme_serif

# ------------------------------------------------------------------------------
# Panel B: Dung removal by D. gazella
# ------------------------------------------------------------------------------
plotRe <- ggplot(GLMMDB, aes(x = BurningRegime, y = RemovalSpecies,
                             fill = BurningRegime)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.6) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1),
              size = 2, alpha = 0.4, color = "cyan4") +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  scale_fill_manual(values = fill_vals) +
  ylab(expression(bold(atop("Dung removal by",
                            paste(bolditalic("D. gazella"), " (g)"))))) +
  xlab("Burning regime") +
  # Post-hoc letter annotations
  annotate("text", label = "A", size = 4, x = 1,   y =   5, family = "serif") +
  annotate("text", label = "A", size = 4, x = 2,   y = 100, family = "serif") +
  annotate("text", label = "A", size = 4, x = 3,   y =  95, family = "serif") +
  annotate("text", label = "B", size = 6, x = 3.5, y = 114, family = "serif") +
  theme_serif

# Combine panels and export
ggarrange(plotReT, plotRe, ncol = 2, nrow = 1)
ggsave("figures/RemovalDung.tiff", plotReT,
       width = 13, height = 9, dpi = 600, bg = "white")

# ==============================================================================
# FIGURE 3: Functional group abundances
# ==============================================================================

# ------------------------------------------------------------------------------
# Panel A: Food relocation strategy (Roller vs. Tunneler)
# ------------------------------------------------------------------------------
FoodR <- GLMMDB |>
  pivot_longer(cols = c(Tunnelers, Rollers),
               names_to  = "FoodRelocation",
               values_to = "Total") |>
  mutate(FoodRelocation = factor(FoodRelocation,
                                 levels = c("Rollers", "Tunnelers"),
                                 labels = c("Roller", "Tunneler"))) |>
  ggplot(aes(x = FoodRelocation, y = Total, fill = BurningRegime)) +
  geom_boxplot(width = 0.6, alpha = 0.6) +
  scale_fill_manual(values = fill_vals) +
  scale_y_continuous(breaks = seq(0, 120, by = 20), limits = c(0, 120)) +
  annotate("text", label = "a", size = 6, x = 0.8, y =  63, family = "serif") +
  annotate("text", label = "a", size = 6, x = 1.0, y =  98, family = "serif") +
  annotate("text", label = "a", size = 6, x = 1.2, y =  83, family = "serif") +
  annotate("text", label = "a", size = 6, x = 1.8, y =  55, family = "serif") +
  annotate("text", label = "b", size = 6, x = 2.0, y =  33, family = "serif") +
  annotate("text", label = "a", size = 6, x = 2.2, y =  83, family = "serif") +
  labs(x     = "Food relocation strategy",
       y     = "Number of individuals",
       fill  = "Burning regime",
       color = "Burning regime") +
  theme_serif2

# ------------------------------------------------------------------------------
# Panel B: Body size (Small, Medium, Large)
# ------------------------------------------------------------------------------
Size <- GLMMDB |>
  pivot_longer(cols = c(Small, Medium, Large),
               names_to  = "Size",
               values_to = "Total") |>
  mutate(Size = factor(Size, levels = c("Small", "Medium", "Large"))) |>
  ggplot(aes(x = Size, y = Total, fill = BurningRegime)) +
  geom_boxplot2(width = 0.6, width.errorbar = 0.1, alpha = 0.6) +
  scale_fill_manual(values = fill_vals) +
  scale_y_continuous(breaks = seq(0, 120, by = 20), limits = c(0, 120)) +
  annotate("text", label = "a", size = 6, x = 0.8, y =  60, family = "serif") +
  annotate("text", label = "a", size = 6, x = 1.0, y =  23, family = "serif") +
  annotate("text", label = "a", size = 6, x = 1.2, y =  70, family = "serif") +
  annotate("text", label = "a", size = 6, x = 1.8, y =  60, family = "serif") +
  annotate("text", label = "a", size = 6, x = 2.0, y =  53, family = "serif") +
  annotate("text", label = "a", size = 6, x = 2.2, y =  65, family = "serif") +
  annotate("text", label = "a", size = 6, x = 2.8, y =  23, family = "serif") +
  annotate("text", label = "a", size = 6, x = 3.0, y =  21, family = "serif") +
  annotate("text", label = "a", size = 6, x = 3.2, y =  20, family = "serif") +
  labs(x     = "Size",
       y     = "Number of individuals",
       fill  = "Burning regime",
       color = "Burning regime") +
  theme_serif2

# ------------------------------------------------------------------------------
# Panel C: Feeding strategy (Coprophagous, Necrophagous, Generalist)
# ------------------------------------------------------------------------------
FS <- GLMMDB |>
  pivot_longer(cols = c(FeedStratC, FeedStratN, FeedStratG),
               names_to  = "FeedingStrategy",
               values_to = "Total") |>
  mutate(FeedingStrategy = factor(
    FeedingStrategy,
    levels = c("FeedStratC", "FeedStratN", "FeedStratG"),
    labels = c("Coprophagous", "Necrophagous", "Generalist")
  )) |>
  ggplot(aes(x = FeedingStrategy, y = Total, fill = BurningRegime)) +
  geom_boxplot2(width = 0.6, width.errorbar = 0.1, alpha = 0.6) +
  scale_fill_manual(values = fill_vals) +
  scale_y_continuous(breaks = seq(0, 120, by = 20), limits = c(0, 125)) +
  annotate("text", label = "ab", size = 6, x = 0.8, y =  53, family = "serif") +
  annotate("text", label = "a",  size = 6, x = 1.0, y =  35, family = "serif") +
  annotate("text", label = "b",  size = 6, x = 1.2, y = 118, family = "serif") +
  annotate("text", label = "a",  size = 6, x = 1.8, y =  55, family = "serif") +
  annotate("text", label = "ab", size = 6, x = 2.0, y =  33, family = "serif") +
  annotate("text", label = "b",  size = 6, x = 2.2, y =  25, family = "serif") +
  annotate("text", label = "a",  size = 6, x = 2.8, y =  45, family = "serif") +
  annotate("text", label = "a",  size = 6, x = 3.0, y =  25, family = "serif") +
  annotate("text", label = "a",  size = 6, x = 3.2, y =  20, family = "serif") +
  labs(x     = "Feeding strategy",
       y     = "Number of individuals",
       fill  = "Burning regime",
       color = "Burning regime") +
  theme_serif2

# ------------------------------------------------------------------------------
# Panel D: Habitat preference (Specialist, Generalist)
# ------------------------------------------------------------------------------
HP <- GLMMDB |>
  pivot_longer(cols = c(HabPrefS, HabPrefG),
               names_to  = "HabitatPreference",
               values_to = "Total") |>
  mutate(HabitatPreference = factor(
    HabitatPreference,
    levels = c("HabPrefS", "HabPrefG"),
    labels = c("Specialist", "Generalist")
  )) |>
  ggplot(aes(x = HabitatPreference, y = Total, fill = BurningRegime)) +
  geom_boxplot2(width = 0.6, width.errorbar = 0.1, alpha = 0.6) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0),
              size = 2, alpha = 0.4, color = "cyan4") +
  scale_fill_manual(values = fill_vals) +
  scale_y_continuous(breaks = seq(0, 120, by = 20), limits = c(0, 125)) +
  annotate("text", label = "a", size = 6, x = 0.8, y =  45, family = "serif") +
  annotate("text", label = "a", size = 6, x = 1.0, y =  30, family = "serif") +
  annotate("text", label = "a", size = 6, x = 1.2, y = 125, family = "serif") +
  annotate("text", label = "a", size = 6, x = 1.8, y = 125, family = "serif") +
  annotate("text", label = "a", size = 6, x = 2.0, y =  35, family = "serif") +
  annotate("text", label = "a", size = 6, x = 2.2, y =  25, family = "serif") +
  labs(x     = "Habitat preference",
       y     = "Number of individuals",
       fill  = "Burning regime",
       color = "Burning regime") +
  theme_serif2 +
  theme(plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm"))

# ------------------------------------------------------------------------------
# Assemble and export Figure 3
# ------------------------------------------------------------------------------
fiGroups <- ggarrange(
  FoodR,
  Size  + theme(legend.position = "none"),
  FS    + theme(legend.position = "none"),
  HP    + theme(legend.position = "none"),
  ncol = 2, nrow = 2,
  align        = "hv",
  common.legend = TRUE,
  legend       = "top"
)

annotate_figure(
  fiGroups,
  left = text_grob("Number of individuals", rot = 90, vjust = 0.7,
                   size = 20, face = "bold", family = "serif")
)

ggsave("figures/FunctionalGroups.svg", width = 13, height = 9,
       dpi = 600, bg = "white")

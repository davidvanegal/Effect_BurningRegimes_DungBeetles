# ==============================================================================
# FuncDiversity.R
# Functional composition (CWM) and proportional abundance by trait
# Study: Dung beetle functional diversity across burning regimes
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load required libraries
# ------------------------------------------------------------------------------
library(tidyverse)      # Data manipulation and visualization
library(lme4)           # Mixed-effects models (loaded for downstream analyses)
library(FD)             # Functional diversity indices
library(multcomp)       # Multiple comparisons
library(performance)    # Model diagnostics
library(RVAideMemoire)  # Permutation tests and model checks
library(ggpubr)         # Figure assembly

# ------------------------------------------------------------------------------
# 2. Compute functional composition (CWM)
# ------------------------------------------------------------------------------

# Load trait matrix (species x traits)
GFDB2 <- read.csv("CWM/GroupFunc.csv", sep = ";")[, -1]

# Load species abundance matrix (sites x species)
DB2    <- read.csv("CWM/Dom.csv", sep = ";")
DomMat <- as.matrix(DB2[, -1])

# Assign species names as row names in trait matrix
row.names(GFDB2) <- colnames(DB2[, -1])

# Compute community-weighted means for all categorical traits
DomGF <- functcomp(GFDB2, DomMat, CWM.type = "all")
write.csv(DomGF, "CWM/RelAbund.csv")

# ------------------------------------------------------------------------------
# 3. Define shared palette and theme
# ------------------------------------------------------------------------------

# Color palette per burning regime
colors <- c(
  "Null" = "lightgoldenrod1",
  "Low"  = "orange",
  "High" = "orangered"
)

# Shared plot theme
theme_serif <- theme_classic(base_family = "serif") +
  theme(
    axis.title.x  = element_text(size = 20, margin = margin(t = 10), face = "bold"),
    axis.title.y  = element_blank(),
    axis.text.x   = element_text(size = 18),
    axis.text.y   = element_text(size = 18),
    legend.title  = element_text(size = 20, face = "bold"),
    legend.text   = element_text(size = 18),
    panel.border  = element_blank(),
    panel.grid    = element_blank(),
    panel.background = element_blank(),
    axis.line     = element_line(color = "black")
  )

# ------------------------------------------------------------------------------
# 4. Generic plotting function for proportional trait abundances
# ------------------------------------------------------------------------------

#' plot_trait
#'
#' Produces a stacked bar chart of proportional abundances for a set of traits.
#'
#' @param data   A data.frame containing trait columns and BurningRegime
#' @param traits Character vector of column names to include
#' @param label  X-axis label for the plot
#' @return A ggplot object

plot_trait <- function(data, traits, label) {
  data |>
    pivot_longer(cols = traits, names_to = "Trait", values_to = "Total") |>
    ggplot(aes(x = Trait, y = Total, fill = BurningRegime)) +
    geom_col(position = "fill", width = 0.6) +
    scale_fill_manual(values = colors) +
    labs(x = label, fill = "Burning regime") +
    theme_serif +
    theme(legend.position = "none")
}

# ------------------------------------------------------------------------------
# 5. Load and prepare relative abundance data
# ------------------------------------------------------------------------------
RelDB <- read.csv("CWM/RelAbund.csv", sep = ";") |>
  mutate(BurningRegime = factor(BurningRegime, levels = c("Null", "Low", "High")))

# Rename trait columns to interpretable labels
RelDB <- RelDB |>
  rename(
    Tunneler          = RelocationStrategy_Tunellers,
    Roller            = RelocationStrategy_Rollers,
    Small             = Size_SmallSize,
    Medium            = Size_MediumSize,
    Large             = Size_LargeSize,
    Coprophagous      = FeedingStrategy_FeedingStrategyC,
    Necrophagous      = FeedingStrategy_FeedingStrategyN,
    Generalist        = FeedingStrategy_FeedingStrategyG,
    GeneralistHabitat = HabitatPrefence_HabitatPrefenceG,
    Specialist        = HabitatPrefence_HabitatPrefenceS
  )

# ------------------------------------------------------------------------------
# 6. Build individual trait panels
# ------------------------------------------------------------------------------

# Food relocation strategy
FRPlot <- plot_trait(RelDB, c("Tunneler", "Roller"),
                     label = "Food relocation strategy")

# Body size
SPlot <- plot_trait(RelDB, c("Small", "Medium", "Large"),
                    label = "Size")

# Feeding strategy (with significance annotations)
FSPlot <- plot_trait(RelDB, c("Coprophagous", "Necrophagous", "Generalist"),
                     label = "Feeding strategy") +
  annotate("text", label = "*",  size = 5, x = 1, y = 0.25, family = "serif") +
  annotate("text", label = "\u00B0", size = 5, x = 1, y = 0.87, family = "serif") +
  annotate("text", label = "\u00B0", size = 5, x = 3, y = 0.06, family = "serif") +
  annotate("text", label = "*",  size = 5, x = 3, y = 0.70, family = "serif")

# Habitat preference (rename Generalist column to avoid conflict)
RelDB2 <- RelDB |>
  dplyr::select(-Generalist) |>
  dplyr::rename(Generalist = GeneralistHabitat)

HPPlot <- plot_trait(RelDB2, c("Generalist", "Specialist"),
                     label = "Habitat preference") +
  annotate("text", label = "*",  x = 1, y = 0.72, size = 5, family = "serif") +
  annotate("text", label = "\u00B0", x = 1, y = 0.07, size = 5, family = "serif") +
  annotate("text", label = "\u00B0", x = 2, y = 0.90, size = 5, family = "serif") +
  annotate("text", label = "*",  x = 2, y = 0.22, size = 5, family = "serif")

# ------------------------------------------------------------------------------
# 7. Assemble and export final figure
# ------------------------------------------------------------------------------
fiGroups <- ggarrange(
  FRPlot,
  SPlot,
  FSPlot + theme(legend.position = "top"),
  HPPlot,
  ncol = 2, nrow = 2,
  align        = "hv",
  common.legend = TRUE,
  legend       = "top"
)

annotate_figure(
  fiGroups,
  left = text_grob("Proportional abundance", rot = 90, vjust = 0.7,
                   size = 20, face = "bold", family = "serif")
)

ggsave("figures/CWM.tiff", width = 11, height = 7, dpi = 600, bg = "white")

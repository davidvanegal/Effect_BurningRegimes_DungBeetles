# ==============================================================================
# NMDS.R
# Non-metric Multidimensional Scaling (NMDS) and community composition tests
# Study: Dung beetle assemblages across burning regimes
# ==============================================================================
# Description:
#   Performs NMDS ordination of dung beetle assemblages across three burning
#   regimes (Null, Low, High). Includes PERMANOVA (global and pairwise),
#   betadisper homogeneity of dispersion test, and Tukey comparisons.
#
# Input:
#   NMDS/NMDSBurning.csv — site x species abundance matrix with metadata
#
# Output:
#   NMDS/NMDS.tiff — publication-quality NMDS ordination plot (600 dpi)
#   Console output — PERMANOVA tables, betadisper results, Tukey comparisons
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load required libraries
# ------------------------------------------------------------------------------
libs <- c("vegan", "ggplot2", "dplyr", "ggrepel", "ggforce",
          "multcomp", "pairwiseAdonis", "indicspecies")
invisible(lapply(libs, library, character.only = TRUE))

# ------------------------------------------------------------------------------
# 2. Load and prepare data
# ------------------------------------------------------------------------------
DBData       <- read.csv("NMDS/NMDSBurning.csv", sep = ";")
DB           <- DBData[, 4:21]             # Species abundance columns
BurningRegime <- DBData$BurningRegime      # Treatment factor

# ------------------------------------------------------------------------------
# 3. NMDS ordination (Bray-Curtis dissimilarity, k = 2)
# ------------------------------------------------------------------------------
nmdsDB     <- metaMDS(DB, k = 2, distance = "bray")

# Extract site scores
nmdsdataDB       <- scores(nmdsDB)$sites %>% as.data.frame()
nmdsdataDB$site  <- DBData$Sites
nmdsdataDB$group <- factor(
  rep(c("High", "Low", "Null"), each = 5),
  levels = c("Null", "Low", "High")
)

# ------------------------------------------------------------------------------
# 4. Compute group centroids and spoke segments
# ------------------------------------------------------------------------------
centroidesDB <- nmdsdataDB |>
  group_by(group) |>
  summarise(across(c(NMDS1, NMDS2), mean), .groups = "drop")

# Join point coordinates with centroid coordinates for spoke lines
segmentosDB <- left_join(nmdsdataDB, centroidesDB,
                         by = "group", suffix = c(".x", ".y"))

# ------------------------------------------------------------------------------
# 5. Define custom plot theme
# ------------------------------------------------------------------------------
theme_custom <- theme_minimal(base_family = "serif") +
  theme(
    axis.title       = element_text(size = 20, face = "bold"),
    axis.text        = element_text(size = 18),
    legend.text      = element_text(size = 18),
    legend.title     = element_text(size = 20, face = "bold"),
    legend.position  = "top",
    axis.line        = element_line(colour = "black"),
    panel.grid       = element_blank(),
    axis.title.x     = element_text(margin = margin(t = 10)),
    axis.title.y     = element_text(margin = margin(t = 10)),
    panel.grid.major = element_line(color = "gray90", linetype = "dashed")
  )

# Color and shape palettes by burning regime
fill_palette  <- c("Null" = "lightgoldenrod2", "Low" = "orange",  "High" = "orangered")
color_palette <- c("Null" = "lightgoldenrod2", "Low" = "orange",  "High" = "orangered")
shape_palette <- c("Null" = 21,                "Low" = 22,        "High" = 24)

# ------------------------------------------------------------------------------
# 6. NMDS ordination plot
# ------------------------------------------------------------------------------
ggplot(nmdsdataDB, aes(x = NMDS1, y = NMDS2)) +

  # Spoke lines from points to group centroids
  geom_segment(
    data    = segmentosDB,
    aes(x = NMDS1.x, y = NMDS2.x, xend = NMDS1.y, yend = NMDS2.y,
        color = group),
    linetype  = "dashed",
    linewidth = 0.5
  ) +

  # Individual site points
  geom_point(aes(shape = group, fill = group),
             size = 5, color = "black", stroke = 1) +

  # Group centroid points (hollow)
  geom_point(data  = centroidesDB,
             aes(x = NMDS1, y = NMDS2, shape = group),
             size = 5, color = "black", stroke = 1.5, fill = "white") +

  # Confidence ellipses per group
  geom_mark_ellipse(expand = 0, aes(color = group, fill = group)) +

  # Site labels (repelled to avoid overlap)
  geom_label_repel(aes(label = site), family = "serif", size = 5) +

  ylim(-0.75, 0.75) +
  scale_fill_manual(values  = fill_palette) +
  scale_color_manual(values = color_palette) +
  scale_shape_manual(values = shape_palette) +
  labs(
    x     = "NMDS1",
    y     = "NMDS2",
    color = "Burning regime",
    shape = "Burning regime",
    fill  = "Burning regime"
  ) +
  theme_custom

ggsave("figures/NMDS.tiff", width = 13, height = 9, dpi = 600, bg = "white")

# ------------------------------------------------------------------------------
# 7. PERMANOVA — global and pairwise
# ------------------------------------------------------------------------------

# Hellinger transformation then Bray-Curtis dissimilarity
DB.mat  <- sqrt(as.matrix(DB))
DB.dist <- vegdist(DB.mat, method = "bray")

# Global PERMANOVA
DB.div <- adonis2(DB.dist ~ as.factor(BurningRegime),
                  data = DBData, permutations = 9999)
print(DB.div)

# Pairwise PERMANOVA
pair_perm <- pairwise.adonis2(DB.dist ~ BurningRegime,
                              data = DBData, permutations = 9999)
print(pair_perm)

# ------------------------------------------------------------------------------
# 8. Betadisper — homogeneity of multivariate dispersion
# ------------------------------------------------------------------------------
betadisper_result <- betadisper(DB.dist, as.factor(DBData$BurningRegime))

# Permutation test for dispersion differences
permutest(betadisper_result)

# Visual inspection of dispersion
plot(betadisper_result)  # Standard deviation ellipses

# Tukey HSD comparisons of mean distances to centroid
tukey_result <- TukeyHSD(betadisper_result)
print(tukey_result)

# Boxplot of distances to centroid by burning regime
boxplot(betadisper_result$distances ~ DBData$BurningRegime)

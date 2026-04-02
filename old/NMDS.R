# 📦 Load essential packages
libs <- c("vegan", "ggplot2", "dplyr", "ggrepel", "ggforce",
          "multcomp", "pairwiseAdonis", "indicspecies")
invisible(lapply(libs, library, character.only = TRUE))

# 📂 Read and prepare data
DBData <- read.csv("NMDS/NMDSBurning.csv", sep = ";")
DB <- DBData[, 4:21]
BurningRegime <- DBData$BurningRegime

# 🔄 NMDS analysis
nmdsDB <- metaMDS(DB, k = 2, distance = "bray")
nmdsdataDB <- scores(nmdsDB)$sites %>% as.data.frame()
nmdsdataDB$site <- DBData$Sites
nmdsdataDB$group <- factor(rep(c("High", "Low", "Null"), each = 5), 
                           levels = c("Null", "Low", "High"))

# 📍 Centroids calculation
centroidesDB <- nmdsdataDB |> 
  group_by(group) |> 
  summarise(across(c(NMDS1, NMDS2), mean), .groups = "drop")

# 🔗 Create segments from each point to its group centroid
segmentosDB <- left_join(nmdsdataDB, centroidesDB, by = "group", suffix = c(".x", ".y"))

# 🎨 Custom theme
theme_custom <- theme_minimal(base_family = "serif") +
  theme(axis.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20, face = "bold"),
        legend.position = "top",
        axis.line = element_line(colour = "black"),
        panel.grid = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(t = 10)),
        panel.grid.major = element_line(color = "gray90", 
                                        linetype = "dashed"))

# 📊 NMDS Plot
ggplot(nmdsdataDB, aes(x = NMDS1, y = NMDS2)) +
  geom_segment(data = segmentosDB, aes(x = NMDS1.x, y = NMDS2.x,
                                       xend = NMDS1.y, yend = NMDS2.y,
                                       color = group), linetype = "dashed", 
               linewidth = 0.5) +
  geom_point(aes(shape = group, fill = group), size = 5, color = "black", stroke = 1) +
  geom_point(data = centroidesDB, aes(x = NMDS1, y = NMDS2, shape = group), 
             size = 5, color = "black", stroke = 1.5, fill = "white") +
  geom_mark_ellipse(expand = 0, aes(color = group, fill = group)) +
  
  geom_label_repel(aes(label = site), family = "serif", size = 5) +
  ylim(-0.75, .75)+
  scale_fill_manual(values = c("Null" = "lightgoldenrod2", 
                               "Low"  = "orange", "High" = "orangered")) +
  scale_color_manual(values = c("Null" = "lightgoldenrod2", 
                                "Low"  = "orange", "High" = "orangered")) +
  scale_shape_manual(values = c("Null" = 21, "Low" = 22, "High" = 24)) +
  labs(x = "NMDS1", y = "NMDS2",
       color = "Burning regime", shape = "Burning regime", fill = "Burning regime") +
  theme_custom

ggsave("NMDS/NMDS.tiff",
       width = 13, height = 9, dpi = 600, bg = "white")

# 🧪 PERMANOVA
DB.mat <- sqrt(as.matrix(DB))
DB.dist <- vegdist(DB.mat, method = "bray")

# Global PERMANOVA
DB.div <- adonis2(DB.dist ~ as.factor(BurningRegime), data = DBData, permutations = 9999)

# Pairwise PERMANOVA
pair_perm <- pairwise.adonis2(DB.dist ~ BurningRegime, data = DBData, permutations = 9999)

# 🎯 Betadisper
betadisper_result <- betadisper(DB.dist, as.factor(DBData$BurningRegime))
permutest(betadisper_result)
plot(betadisper_result)  # standard deviation ellipse

# 🔬 Tukey test for dispersion differences
tukey_result <- TukeyHSD(betadisper_result)
boxplot(betadisper_result$distances ~ DBData$BurningRegime)
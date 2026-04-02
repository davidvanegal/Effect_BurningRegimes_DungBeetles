# 📦 Libraries
library(tidyverse)
library(lme4)
library(FD)
library(multcomp)
library(performance)
library(RVAideMemoire)
library(ggpubr)

# 📊 Functional composition

# 📁 Data Import
GFDB2 <- read.csv("CWM/GroupFunc.csv", sep = ";")[,-1]
DB2 <- read.csv("CWM/Dom.csv", sep = ";")

DomMat <- as.matrix(DB2[,-1])
row.names(GFDB2) <- colnames(DB2[,-1])
DomGF <- functcomp(GFDB2, DomMat, CWM.type = "all")
write.csv(DomGF, "CWM/RelAbund.csv")

# 📈 Plot wrapper
# 🎨 Palette & Theme
colors <- c("Null" = "lightgoldenrod1", "Low" = "orange", "High" = "orangered")
theme_serif <- theme_classic(base_family = "serif") +
  theme(axis.title.x = element_text(size = 20, margin = margin(t = 10), face = "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 18),
        panel.border = element_blank(), panel.grid = element_blank(),
        panel.background = element_blank(), axis.line = element_line(color = "black"))

plot_trait <- function(data, traits, label) {
  data |>
    pivot_longer(cols = traits, names_to = "Trait", values_to = "Total") |>
    ggplot(aes(x = Trait, y = Total, fill = BurningRegime)) +
    geom_col(position = "fill", width = 0.6) +
    scale_fill_manual(values = colors) +
    labs(x = label, fill = "Burning regime") +
    theme_serif + theme(legend.position = "none")
}

# 📁 Load RelAbund
RelDB <- read.csv("CWM/RelAbund.csv", sep = ";") |>
  mutate(BurningRegime = factor(BurningRegime, levels = c("Null", "Low", "High")))

RelDB <- RelDB |>
  rename(Tunneler = RelocationStrategy_Tunellers,
    Roller = RelocationStrategy_Rollers,
    Small = Size_SmallSize,
    Medium = Size_MediumSize,
    Large = Size_LargeSize,
    Coprophagous = FeedingStrategy_FeedingStrategyC,
    Necrophagous = FeedingStrategy_FeedingStrategyN,
    Generalist = FeedingStrategy_FeedingStrategyG,
    GeneralistHabitat = HabitatPrefence_HabitatPrefenceG,
    Specialist = HabitatPrefence_HabitatPrefenceS)

# Individual plots
FRPlot <- plot_trait(RelDB, c("Tunneler", "Roller"), "Food relocation strategy")
SPlot  <- plot_trait(RelDB, c("Small", "Medium", "Large"), "Size")
FSPlot <- plot_trait(RelDB, c("Coprophagous", "Necrophagous", 
                              "Generalist"), "Feeding strategy")+
  annotate("text", label = "*", size = 5, x = 1, y = 0.25, family = "serif")+ 
  annotate("text", label = "°", size = 5, x = 1, y = 0.87, family = "serif")+ 
  annotate("text", label = "°", size = 5, x = 3, y = 0.06, family = "serif")+ 
  annotate("text", label = "*", size = 5, x = 3, y = 0.7, family = "serif")

RelDB2 <- RelDB |>
  dplyr::select(-Generalist) |>
  dplyr::rename(Generalist = GeneralistHabitat)


HPPlot <- plot_trait(RelDB2, c("Generalist", "Specialist"), "Habitat preference")+
  annotate("text", label = "*", x = 1, y = 0.72, size = 5, family = "serif") +
  annotate("text", label = "°", x = 1, y = 0.07, size = 5, family = "serif") +
  annotate("text", label = "°", x = 2, y = 0.9, size = 5, family = "serif") +
  annotate("text", label = "*", x = 2, y = 0.22, size = 5, family = "serif")
  

# Final figure
fiGroups <- ggarrange(FRPlot, SPlot, FSPlot + theme(legend.position = "top"), HPPlot,
                      ncol = 2, nrow = 2, align = "hv", common.legend = TRUE, legend = "top")

annotate_figure(fiGroups,
                left = text_grob("Proportional abundance", rot = 90, vjust = 0.7,
                                 size = 20, face = "bold", family = "serif"))

# Save
ggsave("CWM/CWM.tiff", width = 11, height = 7, dpi = 600, bg = "white")

source("RangeAbundance/SpecDist.R")
source("RangeAbundance/DetAbu.R") 
source("RangeAbundance/UndAbu.R") 

## Load packages ----
library(iNEXT)
library(Rcpp)
library(ggplot2)
library(ggrepel)
library(stringr)
library(dplyr)
library(readxl)
library(ggpubr)

## Read data ----
ra <- read.csv("RangeAbundance/RangeAbundance.csv", sep = ";")

ra2 <- ra |> 
  group_by(Species) |> 
  summarise(across(c(Null, Low, High), sum))

## PA ---- 
out1 <- SpecDist(ra2$Null, "abundance") 
out1 <- subset(out1, probability>0) 
totalNull <- dim(out1)[1] 
i=1 
for (i in 1:totalNull) { 
  out1$number[i] <- i 
  out1$sp[i] <- as.numeric(rownames(out1))[i] } 
out1$BurningRegimen <- "Null" 

## 
out2 <- SpecDist(ra2$Low, "abundance") 
out2 <- subset(out2, probability>0) 
totalLow <- dim(out2)[1] 
i=1 
for (i in 1:totalLow) { 
  out2$number[i] <- i 
  out2$sp[i] <- as.numeric(rownames(out2))[i] } 
out2$BurningRegimen <- "Low" 

## VS ---- 
out3 <- SpecDist(ra2$High, "abundance") 
out3 <- subset(out3, probability>0) 
totalHigh <- dim(out3)[1] 
i=1 
for (i in 1:totalHigh) { 
  out3$number[i] <- i 
  out3$sp[i] <- as.numeric(rownames(out3))[i] } 
out3$BurningRegimen <- "High" 

# NULL
out1 <- SpecDist(ra2$Null, "abundance") %>% subset(probability > 0)
out1$number <- seq_len(nrow(out1))
out1$sp <- as.numeric(rownames(out1))
out1$Species <- ra2$Species[out1$sp]
out1$BurningRegimen <- "Null"

# LOW
out2 <- SpecDist(ra2$Low, "abundance") %>% subset(probability > 0)
out2$number <- seq_len(nrow(out2))
out2$sp <- as.numeric(rownames(out2))
out2$Species <- ra2$Species[out2$sp]
out2$BurningRegimen <- "Low"

# HIGH
out3 <- SpecDist(ra2$High, "abundance") %>% subset(probability > 0)
out3$number <- seq_len(nrow(out3))
out3$sp <- as.numeric(rownames(out3))
out3$Species <- ra2$Species[out3$sp]
out3$BurningRegimen <- "High"


## Create a single base 
outtotal <- rbind(out1, out2, out3) 
outtotal <-subset(outtotal, method %in% c("detected")) 
outtotal$BurningRegimen <- factor(outtotal$BurningRegimen, 
                                  levels = c("Null", "Low", "High"))

theme_custom <- theme_minimal(base_family = "serif") +
  theme(axis.title = element_text(size = 20, face = "bold"),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20, face = "bold"),
        legend.position = "top",
        axis.line = element_line(colour = "black"),
        panel.grid = element_blank(),
        axis.title.x = element_text(margin = margin(t = 20)),
        axis.title.y = element_blank(),
        panel.grid.major = element_line(color = "gray95", 
                                        linetype = "dashed"))

## Plot ----
p <- ggplot(outtotal, aes(x = number, y = probability, color = BurningRegimen, 
                          shape = BurningRegimen)) + 
  geom_point(size = 4) + 
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("lightgoldenrod2", "orange", "orangered")) +
  scale_x_continuous("Species rank") +
  scale_y_continuous("Proportional abundance") +
  facet_grid(cols = vars(BurningRegimen)) +
  geom_text(aes(label = Species), size = 5, family = "serif", color = "black", 
            angle = 45, hjust = -0.15, vjust = -0.1, fontface = "bold.italic")+
  ylim(0, .6)+
  xlim(0, 19)+
  labs(col = "Burning regime", shape = "Burning regime")+
  xlab("Species rank")+
  theme_custom

p <- annotate_figure(p, 
                left = text_grob("Number of individuals", rot = 90, vjust = 0.3, 
                                 size = 20, face = "bold", family = "serif"))

## Save plot ----
ggsave("RangeAbundance/RangeAbundance.tiff",
       width = 15, height = 9, dpi = 600, bg = "white")

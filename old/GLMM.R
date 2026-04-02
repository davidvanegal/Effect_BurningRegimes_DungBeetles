library(multcomp) # Simultaneous Inference in General Parametric Models CRAN v1.4-25
library(lme4) # Linear Mixed-Effects Models using 'Eigen' and S4 CRAN v1.1-36 # Linear Mixed-Effects Models using 'Eigen' and S4 CRAN v1.1-35.2
library(lmtest) # Testing Linear Regression Models CRAN v0.9-40
library(lattice) # Trellis Graphics for R CRAN v0.21-9 
library(tidyverse) # Easily Install and Load the 'Tidyverse' CRAN v2.0.0 
library(emmeans) # Estimated Marginal Means, aka Least-Squares Means CRAN v1.10.1
library(psych) # Procedures for Psychological, Psychometric, and Personality Research CRAN v2.4.3
library(performance) # Assessment of Regression Models Performance CRAN v0.13.0
library(MASS) # Support Functions and Datasets for Venables and Ripley's MASS CRAN v7.3-60
library(lme4) # Linear Mixed-Effects Models using 'Eigen' and S4 CRAN v1.1-36
library(glmmTMB)
library(ggeffects)
library(ggplot2)

# Read data
GLMMDB <- read.csv("GLMM/Data.csv",  sep = ";") 
head(GLMMDB)

# Convert factor -----
GLMMDB$BurningRegime <- factor(GLMMDB$BurningRegime, levels = c("Null", 
                                                                "Low",
                                                                "High"))
GLMMDB$GrazingEnviroments <- factor(GLMMDB$GrazingEnviroments, 
                                    levels = c("Secondary forest", 
                                               "Biodiverse pastures",
                                               "Grass monoculture"))

theme_serif <- theme_classic(base_family = "serif") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        panel.grid.major = element_line(color = "gray95", linetype = "dashed"),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        panel.background = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18, face = "bold"))

# Abundance ----
m_abund <- glmmTMB(Abundance ~ BurningRegime + GrazingEnviroments + (1 | IDT),
  family = nbinom2, data = GLMMDB)
summary(m_abund)

## Inflación de ceros
simres_m_abund <- simulateResiduals(m_abund)
testZeroInflation(simres_m_abund)  # Test inflación de ceros

# Predicts
preds <- ggpredict(m_abund, terms = c("BurningRegime", "GrazingEnviroments"),
                   bias_correction = TRUE)
preds$group <- factor(preds$group,
                      levels = c("Secondary forest",
                                 "Biodiverse pastures",
                                 "Grass monoculture"))

preds <- subset(preds, !(x == "High" & group == "Secondary forest"))


# Gráfico combinado
plot_Abund <- ggplot(GLMMDB, aes(x = BurningRegime, y = Abundance, fill = GrazingEnviroments)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.8), colour = "cyan4")+ 
  
  geom_point(data = preds, aes(x = x, y = predicted, color = group),
             position = position_dodge(width = 0.8), size = 3, inherit.aes = FALSE) +
  
  geom_errorbar(data = preds,
                aes(x = x, ymin = conf.low, ymax = conf.high, color = group),
                position = position_dodge(width = 0.8),
                width = 0.5,        # ancho de las "patitas" horizontales
                size = 0.6,         # grosor de la línea vertical
                inherit.aes = FALSE)+
  
  scale_fill_manual(values = c("lightgoldenrod1", "orange", "orangered")) +
  scale_color_manual(values = c("lightgoldenrod1", "orange", "orangered")) +
  
  labs(y = "Abundance (Number of individuals)", x = "Burning Regime",
       fill = "Grazing environments", color = "Grazing environments") +
  theme_serif

# Roller ----
m_roller <- glmmTMB(Roller ~ BurningRegime + GrazingEnviroments + (1 | IDT),
                   family = nbinom2, data = GLMMDB)
summary(m_roller)

## Inflación de ceros
simres_roller <- simulateResiduals(m_roller)
testZeroInflation(simres_roller)  # Test inflación de ceros

# Predicts
preds_roller <- ggpredict(m_roller, terms = c("BurningRegime", "GrazingEnviroments"),
                   bias_correction = TRUE)
preds_roller$group <- factor(preds_roller$group,
                      levels = c("Secondary forest",
                                 "Biodiverse pastures",
                                 "Grass monoculture"))

preds_roller <- subset(preds_roller, !(x == "High" & group == "Secondary forest"))

# Gráfico combinado
plot_Roller <- ggplot(GLMMDB, aes(x = BurningRegime, y = Roller, fill = GrazingEnviroments)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.8), colour = "cyan4")+ 
  
  geom_point(data = preds_roller, aes(x = x, y = predicted, color = group),
             position = position_dodge(width = 0.8), size = 3, inherit.aes = FALSE) +
  
  geom_errorbar(data = preds_roller,
                aes(x = x, ymin = conf.low, ymax = conf.high, color = group),
                position = position_dodge(width = 0.8),
                width = 0.5,        # ancho de las "patitas" horizontales
                size = 0.6,         # grosor de la línea vertical
                inherit.aes = FALSE)+
  
  scale_fill_manual(values = c("lightgoldenrod1", "orange", "orangered")) +
  scale_color_manual(values = c("lightgoldenrod1", "orange", "orangered")) +
  
  labs(y = "Abundance (Number of individuals)", x = "Burning Regime",
       fill = "Grazing environments", color = "Grazing environments") +
  theme_serif

# Tunneler ----
m_tunneler <- glmmTMB(Tuneller ~ BurningRegime + GrazingEnviroments + (1 | IDT),
                    family = nbinom2, data = GLMMDB)
summary(m_tunneler)

## Inflación de ceros
simres_tunneler <- simulateResiduals(m_tunneler)
testZeroInflation(simres_tunneler)  # Test inflación de ceros

# Predicts
preds_tunneler <- ggpredict(m_tunneler, terms = c("BurningRegime", "GrazingEnviroments"),
                          bias_correction = TRUE)
preds_tunneler$group <- factor(preds_tunneler$group,
                             levels = c("Secondary forest",
                                        "Biodiverse pastures",
                                        "Grass monoculture"))

preds_tunneler <- subset(preds_tunneler, !(x == "High" & group == "Secondary forest"))

# Gráfico combinado
plot_tunneler <- ggplot(GLMMDB, aes(x = BurningRegime, y = Tuneller, fill = GrazingEnviroments)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.8), colour = "cyan4")+ 
  
  geom_point(data = preds_tunneler, aes(x = x, y = predicted, color = group),
             position = position_dodge(width = 0.8), size = 3, inherit.aes = FALSE) +
  
  geom_errorbar(data = preds_tunneler,
                aes(x = x, ymin = conf.low, ymax = conf.high, color = group),
                position = position_dodge(width = 0.8),
                width = 0.5,        # ancho de las "patitas" horizontales
                size = 0.6,         # grosor de la línea vertical
                inherit.aes = FALSE)+
  
  scale_fill_manual(values = c("lightgoldenrod1", "orange", "orangered")) +
  scale_color_manual(values = c("lightgoldenrod1", "orange", "orangered")) +
  
  labs(y = "Abundance (Number of individuals)", x = "Burning Regime",
       fill = "Grazing environments", color = "Grazing environments") +
  theme_serif

# Large ----
m_Large <- glmmTMB(Large ~ BurningRegime + GrazingEnviroments + (1 | IDT),
                      family = nbinom2, data = GLMMDB)
summary(m_Large)

## Inflación de ceros
simres_Large <- simulateResiduals(m_Large)
testZeroInflation(simres_Large)  # Test inflación de ceros

# Predicts
preds_Large <- ggpredict(m_Large, terms = c("BurningRegime", "GrazingEnviroments"),
                            bias_correction = TRUE)
preds_Large$group <- factor(preds_Large$group,
                               levels = c("Secondary forest",
                                          "Biodiverse pastures",
                                          "Grass monoculture"))

preds_Large <- subset(preds_Large, !(x == "High" & group == "Secondary forest"))

# Gráfico combinado
plot_Large <- ggplot(GLMMDB, aes(x = BurningRegime, y = Large, fill = GrazingEnviroments)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.8), colour = "cyan4")+ 
  
  geom_point(data = preds_Large, aes(x = x, y = predicted, color = group),
             position = position_dodge(width = 0.8), size = 3, inherit.aes = FALSE) +
  
  geom_errorbar(data = preds_Large,
                aes(x = x, ymin = conf.low, ymax = conf.high, color = group),
                position = position_dodge(width = 0.8),
                width = 0.5,        # ancho de las "patitas" horizontales
                size = 0.6,         # grosor de la línea vertical
                inherit.aes = FALSE)+
  
  scale_fill_manual(values = c("lightgoldenrod1", "orange", "orangered")) +
  scale_color_manual(values = c("lightgoldenrod1", "orange", "orangered")) +
  
  labs(y = "Abundance (Number of individuals)", x = "Burning Regime",
       fill = "Grazing environments", color = "Grazing environments") +
  theme_serif

# Medium ----
m_Medium <- glmmTMB(Medium ~ BurningRegime + GrazingEnviroments + (1 | IDT),
                   family = nbinom2, data = GLMMDB)
summary(m_Medium)

## Inflación de ceros
simres_Medium <- simulateResiduals(m_Medium)
testZeroInflation(simres_Medium)  # Test inflación de ceros

# Predicts
preds_Medium <- ggpredict(m_Medium, terms = c("BurningRegime", "GrazingEnviroments"),
                         bias_correction = TRUE)
preds_Medium$group <- factor(preds_Medium$group,
                            levels = c("Secondary forest",
                                       "Biodiverse pastures",
                                       "Grass monoculture"))

preds_Medium <- subset(preds_Medium, !(x == "High" & group == "Secondary forest"))

# Gráfico combinado
plot_Medium <- ggplot(GLMMDB, aes(x = BurningRegime, y = Medium, fill = GrazingEnviroments)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.8), colour = "cyan4")+ 
  
  geom_point(data = preds_Medium, aes(x = x, y = predicted, color = group),
             position = position_dodge(width = 0.8), size = 3, inherit.aes = FALSE) +
  
  geom_errorbar(data = preds_Medium,
                aes(x = x, ymin = conf.low, ymax = conf.high, color = group),
                position = position_dodge(width = 0.8),
                width = 0.5,        # ancho de las "patitas" horizontales
                size = 0.6,         # grosor de la línea vertical
                inherit.aes = FALSE)+
  
  scale_fill_manual(values = c("lightgoldenrod1", "orange", "orangered")) +
  scale_color_manual(values = c("lightgoldenrod1", "orange", "orangered")) +
  
  labs(y = "Abundance (Number of individuals)", x = "Burning Regime",
       fill = "Grazing environments", color = "Grazing environments") +
  theme_serif

# Small ----
m_Small <- glmmTMB(Small ~ BurningRegime + GrazingEnviroments + (1 | IDT),
                   family = nbinom2, data = GLMMDB)
summary(m_Small)

## Inflación de ceros
simres_Small <- simulateResiduals(m_Small)
testZeroInflation(simres_Small)  # Test inflación de ceros

# Predicts
preds_Small <- ggpredict(m_Small, terms = c("BurningRegime", "GrazingEnviroments"),
                         bias_correction = TRUE)
preds_Small$group <- factor(preds_Small$group,
                            levels = c("Secondary forest",
                                       "Biodiverse pastures",
                                       "Grass monoculture"))

preds_Small <- subset(preds_Small, !(x == "High" & group == "Secondary forest"))

# Gráfico combinado
plot_Small <- ggplot(GLMMDB, aes(x = BurningRegime, y = Small, fill = GrazingEnviroments)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.8),
    colour = "cyan4", outlier.shape = NA) +
  coord_cartesian(ylim = c(0, quantile(GLMMDB$Small, 0.97)))+
  geom_point(data = preds_Small, aes(x = x, y = predicted, color = group),
             position = position_dodge(width = 0.8), size = 3, inherit.aes = FALSE) +
  
  geom_errorbar(data = preds_Small,
                aes(x = x, ymin = conf.low, ymax = conf.high, color = group),
                position = position_dodge(width = 0.8),
                width = 0.5,        # ancho de las "patitas" horizontales
                size = 0.6,         # grosor de la línea vertical
                inherit.aes = FALSE)+
  
  scale_fill_manual(values = c("lightgoldenrod1", "orange", "orangered")) +
  scale_color_manual(values = c("lightgoldenrod1", "orange", "orangered")) +
  
  labs(y = "Abundance (Number of individuals)", x = "Burning Regime",
       fill = "Grazing environments", color = "Grazing environments") +
  theme_serif

# HabPrefG ----
m_HabPrefG <- glmmTMB(HabPrefG ~ BurningRegime + GrazingEnviroments + (1 | IDT),
                   family = nbinom2, data = GLMMDB)
summary(m_HabPrefG)

## Inflación de ceros
simres_HabPrefG <- simulateResiduals(m_HabPrefG)
testZeroInflation(simres_HabPrefG)  # Test inflación de ceros

# Predicts
preds_HabPrefG <- ggpredict(m_HabPrefG, terms = c("BurningRegime", "GrazingEnviroments"),
                         bias_correction = TRUE)
preds_HabPrefG$group <- factor(preds_HabPrefG$group,
                            levels = c("Secondary forest",
                                       "Biodiverse pastures",
                                       "Grass monoculture"))

preds_HabPrefG <- subset(preds_HabPrefG, !(x == "High" & group == "Secondary forest"))

# Gráfico combinado
plot_HabPrefG <- ggplot(GLMMDB, aes(x = BurningRegime, y = HabPrefG, fill = GrazingEnviroments)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.8),
               colour = "cyan4", outlier.shape = NA) +
  coord_cartesian(ylim = c(0, quantile(GLMMDB$HabPrefG, 0.97)))+
  geom_point(data = preds_HabPrefG, aes(x = x, y = predicted, color = group),
             position = position_dodge(width = 0.8), size = 3, inherit.aes = FALSE) +
  
  geom_errorbar(data = preds_HabPrefG,
                aes(x = x, ymin = conf.low, ymax = conf.high, color = group),
                position = position_dodge(width = 0.8),
                width = 0.5,        # ancho de las "patitas" horizontales
                size = 0.6,         # grosor de la línea vertical
                inherit.aes = FALSE)+
  
  scale_fill_manual(values = c("lightgoldenrod1", "orange", "orangered")) +
  scale_color_manual(values = c("lightgoldenrod1", "orange", "orangered")) +
  
  labs(y = "Abundance (Number of individuals)", x = "Burning Regime",
       fill = "Grazing environments", color = "Grazing environments") +
  theme_serif

# HabPrefS ----
m_HabPrefS <- glmmTMB(HabPrefS ~ BurningRegime + GrazingEnviroments + (1 | IDT),
                      family = nbinom2, data = GLMMDB)
summary(m_HabPrefS)

## Inflación de ceros
simres_HabPrefS <- simulateResiduals(m_HabPrefS)
testZeroInflation(simres_HabPrefS)  # Test inflación de ceros

# Predicts
preds_HabPrefS <- ggpredict(m_HabPrefS, terms = c("BurningRegime", "GrazingEnviroments"),
                            bias_correction = TRUE)
preds_HabPrefS$group <- factor(preds_HabPrefS$group,
                               levels = c("Secondary forest",
                                          "Biodiverse pastures",
                                          "Grass monoculture"))

preds_HabPrefS <- subset(preds_HabPrefS, !(x == "High" & group == "Secondary forest"))

# Gráfico combinado
plot_HabPrefS <- ggplot(GLMMDB, aes(x = BurningRegime, y = HabPrefS, fill = GrazingEnviroments)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.8),
               colour = "cyan4", outlier.shape = NA) +
  coord_cartesian(ylim = c(0, quantile(GLMMDB$HabPrefS, 0.97)))+
  geom_point(data = preds_HabPrefS, aes(x = x, y = predicted, color = group),
             position = position_dodge(width = 0.8), size = 3, inherit.aes = FALSE) +
  
  geom_errorbar(data = preds_HabPrefS,
                aes(x = x, ymin = conf.low, ymax = conf.high, color = group),
                position = position_dodge(width = 0.8),
                width = 0.5,        # ancho de las "patitas" horizontales
                size = 0.6,         # grosor de la línea vertical
                inherit.aes = FALSE)+
  
  scale_fill_manual(values = c("lightgoldenrod1", "orange", "orangered")) +
  scale_color_manual(values = c("lightgoldenrod1", "orange", "orangered")) +
  
  labs(y = "Abundance (Number of individuals)", x = "Burning Regime",
       fill = "Grazing environments", color = "Grazing environments") +
  theme_serif

# FeedC ----
m_FeedC <- glmmTMB(FeedC ~ BurningRegime + GrazingEnviroments + (1 | IDT),
                      family = nbinom2, data = GLMMDB)
summary(m_FeedC)

## Inflación de ceros
simres_FeedC <- simulateResiduals(m_FeedC)
testZeroInflation(simres_FeedC)  # Test inflación de ceros

# Predicts
preds_FeedC <- ggpredict(m_FeedC, terms = c("BurningRegime", "GrazingEnviroments"),
                            bias_correction = TRUE)
preds_FeedC$group <- factor(preds_FeedC$group,
                               levels = c("Secondary forest",
                                          "Biodiverse pastures",
                                          "Grass monoculture"))

preds_FeedC <- subset(preds_FeedC, !(x == "High" & group == "Secondary forest"))

# Gráfico combinado
plot_FeedC <- ggplot(GLMMDB, aes(x = BurningRegime, y = FeedC, fill = GrazingEnviroments)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.8),
               colour = "cyan4", outlier.shape = NA) +
  coord_cartesian(ylim = c(0, quantile(GLMMDB$FeedC, 0.97)))+
  geom_point(data = preds_FeedC, aes(x = x, y = predicted, color = group),
             position = position_dodge(width = 0.8), size = 3, inherit.aes = FALSE) +
  
  geom_errorbar(data = preds_FeedC,
                aes(x = x, ymin = conf.low, ymax = conf.high, color = group),
                position = position_dodge(width = 0.8),
                width = 0.5,        # ancho de las "patitas" horizontales
                size = 0.6,         # grosor de la línea vertical
                inherit.aes = FALSE)+
  
  scale_fill_manual(values = c("lightgoldenrod1", "orange", "orangered")) +
  scale_color_manual(values = c("lightgoldenrod1", "orange", "orangered")) +
  
  labs(y = "Abundance (Number of individuals)", x = "Burning Regime",
       fill = "Grazing environments", color = "Grazing environments") +
  theme_serif

# FeedG ----
m_FeedG <- glmmTMB(FeedG ~ BurningRegime + GrazingEnviroments + (1 | IDT),
                   family = nbinom2, data = GLMMDB)
summary(m_FeedG)

## Inflación de ceros
simres_FeedG <- simulateResiduals(m_FeedG)
testZeroInflation(simres_FeedG)  # Test inflación de ceros

# Predicts
preds_FeedG <- ggpredict(m_FeedG, terms = c("BurningRegime", "GrazingEnviroments"),
                         bias_correction = TRUE)
preds_FeedG$group <- factor(preds_FeedG$group,
                            levels = c("Secondary forest",
                                       "Biodiverse pastures",
                                       "Grass monoculture"))

preds_FeedG <- subset(preds_FeedG, !(x == "High" & group == "Secondary forest"))

# Gráfico combinado
plot_FeedG <- ggplot(GLMMDB, aes(x = BurningRegime, y = FeedG, fill = GrazingEnviroments)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.8),
               colour = "cyan4", outlier.shape = NA) +
  coord_cartesian(ylim = c(0, quantile(GLMMDB$FeedG, 0.97)))+
  geom_point(data = preds_FeedG, aes(x = x, y = predicted, color = group),
             position = position_dodge(width = 0.8), size = 3, inherit.aes = FALSE) +
  
  geom_errorbar(data = preds_FeedG,
                aes(x = x, ymin = conf.low, ymax = conf.high, color = group),
                position = position_dodge(width = 0.8),
                width = 0.5,        # ancho de las "patitas" horizontales
                size = 0.6,         # grosor de la línea vertical
                inherit.aes = FALSE)+
  
  scale_fill_manual(values = c("lightgoldenrod1", "orange", "orangered")) +
  scale_color_manual(values = c("lightgoldenrod1", "orange", "orangered")) +
  
  labs(y = "Abundance (Number of individuals)", x = "Burning Regime",
       fill = "Grazing environments", color = "Grazing environments") +
  theme_serif

# FeedN ----
m_FeedN <- glmmTMB(FeedN ~ BurningRegime + GrazingEnviroments + (1 | IDT),
                   family = nbinom2, data = GLMMDB)
summary(m_FeedN)

## Inflación de ceros
simres_FeedN <- simulateResiduals(m_FeedN)
testZeroInflation(simres_FeedN)  # Test inflación de ceros

# Predicts
preds_FeedN <- ggpredict(m_FeedN, terms = c("BurningRegime", "GrazingEnviroments"),
                         bias_correction = TRUE)
preds_FeedN$group <- factor(preds_FeedN$group,
                            levels = c("Secondary forest",
                                       "Biodiverse pastures",
                                       "Grass monoculture"))

preds_FeedN <- subset(preds_FeedN, !(x == "High" & group == "Secondary forest"))

# Gráfico combinado
plot_FeedN <- ggplot(GLMMDB, aes(x = BurningRegime, y = FeedN, fill = GrazingEnviroments)) +
  geom_boxplot(alpha = 0.5, position = position_dodge(width = 0.8),
               colour = "cyan4", outlier.shape = NA) +
  coord_cartesian(ylim = c(0, quantile(GLMMDB$FeedN, 0.97)))+
  geom_point(data = preds_FeedN, aes(x = x, y = predicted, color = group),
             position = position_dodge(width = 0.8), size = 3, inherit.aes = FALSE) +
  
  geom_errorbar(data = preds_FeedN,
                aes(x = x, ymin = conf.low, ymax = conf.high, color = group),
                position = position_dodge(width = 0.8),
                width = 0.5,        # ancho de las "patitas" horizontales
                size = 0.6,         # grosor de la línea vertical
                inherit.aes = FALSE)+
  
  scale_fill_manual(values = c("lightgoldenrod1", "orange", "orangered")) +
  scale_color_manual(values = c("lightgoldenrod1", "orange", "orangered")) +
  
  labs(y = "Abundance (Number of individuals)", x = "Burning Regime",
       fill = "Grazing environments", color = "Grazing environments") +
  theme_serif

## ____________

# Librerías necesarias
library(ggplot2)
library(dplyr)
library(cowplot)
library(gtable)

# -----------------------------
# 1. Definir orden estricto de variables
# -----------------------------
var_order <- rev(c("Abundance","Roller","Tuneller","Large","Medium","Small",
                   "FeedC","FeedG","FeedN", "HabPrefG","HabPrefS"))

# -----------------------------
# 2. Tablas de efectos
# -----------------------------
burning <- tribble(
  ~Variable,     ~Factor,            ~Effect,
  "Abundance",   "Low", "Negative",
  "Abundance",   "High","Neutral",
  "Roller",      "Low", "Neutral",
  "Roller",      "High","Neutral",
  "Tuneller",    "Low", "Negative",
  "Tuneller",    "High","Neutral",
  "Large",       "Low", "Neutral",
  "Large",       "High","Neutral",
  "Medium",      "Low", "Neutral",
  "Medium",      "High","Neutral",
  "Small",       "Low", "Negative",
  "Small",       "High","Neutral",
  "HabPrefG",    "Low", "Negative",
  "HabPrefG",    "High","Negative",
  "HabPrefS",    "Low", "Neutral",
  "HabPrefS",    "High","Positive",
  "FeedC",       "Low", "Neutral",
  "FeedC",       "High","Tendency",
  "FeedG",       "Low", "Neutral",
  "FeedG",       "High","Neutral",
  "FeedN",       "Low", "Tendency",
  "FeedN",       "High","Negative"
) %>% mutate(Variable = factor(Variable, levels = var_order))

cover <- tribble(
  ~Variable,     ~Factor,                       ~Effect,
  "Abundance",   "Biodiverse pastures",         "Neutral",
  "Abundance",   "Grass monoculture",           "Negative",
  "Roller",      "Biodiverse pastures",         "Neutral",
  "Roller",      "Grass monoculture",           "Negative",
  "Tuneller",    "Biodiverse pastures",         "Neutral",
  "Tuneller",    "Grass monoculture",           "Tendency",
  "Large",       "Biodiverse pastures",         "Tendency",
  "Large",       "Grass monoculture",           "Neutral",
  "Medium",      "Biodiverse pastures",         "Neutral",
  "Medium",      "Grass monoculture",           "Negative",
  "Small",       "Biodiverse pastures",         "Neutral",
  "Small",       "Grass monoculture",           "Negative",
  "HabPrefG",    "Biodiverse pastures",         "Negative",
  "HabPrefG",    "Grass monoculture",           "Negative",
  "HabPrefS",    "Biodiverse pastures",         "Neutral",
  "HabPrefS",    "Grass monoculture",           "Negative",
  "FeedC",       "Biodiverse pastures",         "Neutral",
  "FeedC",       "Grass monoculture",           "Neutral",
  "FeedG",       "Biodiverse pastures",         "Negative",
  "FeedG",       "Grass monoculture",           "Negative",
  "FeedN",       "Biodiverse pastures",         "Neutral",
  "FeedN",       "Grass monoculture",           "Tendency"
) %>% mutate(Variable = factor(Variable, levels = var_order))

# -----------------------------
# 3. Definir colores
# -----------------------------
cols <- c("Negative"="firebrick", "Neutral"="gray95",
          "Positive"="darkgreen", "Tendency"="gold")

# -----------------------------
# 4. Crear gráficos
# -----------------------------
p_burn <- ggplot(burning, aes(x=Factor, y=Variable, fill=Effect)) +
  geom_tile(color="white") +
  scale_fill_manual(values=cols) +
  labs(x="Burning regime", y="Variable") +
  theme_serif

p_cover <- ggplot(cover, aes(x=Factor, y=Variable, fill=Effect)) +
  geom_tile(color="white") +
  scale_fill_manual(values=cols) +
  labs(x="Grazing environment", y=NULL) +
  theme_serif

# -----------------------------
# 5. Función para extraer la leyenda
# -----------------------------
get_legend <- function(myplot) {
  tmp <- ggplotGrob(myplot)
  leg <- gtable::gtable_filter(tmp, "guide-box")
  return(leg)
}

legend <- get_legend(p_burn)

# -----------------------------
# 6. Quitar leyenda de ambos gráficos
# -----------------------------
p_burn_nolegend <- p_burn + theme(legend.position = "none")
p_cover_nolegend <- p_cover + theme(legend.position = "none",
                                    axis.text.y = element_text(color="white"))

# -----------------------------
# 7. Combinar gráficos en horizontal y agregar leyenda única
# -----------------------------
plots <- plot_grid(p_burn_nolegend, p_cover_nolegend, ncol = 2)
combined <- plot_grid(plots, legend, ncol = 1, rel_heights = c(1, 0.15))

# Mostrar resultado
combined

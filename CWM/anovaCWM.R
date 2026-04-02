# ==============================================================================
# anovaCWM.R
# Statistical results — Functional Groups (one-way ANOVA)
# Study: Dung beetle functional trait proportions across burning regimes
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load required libraries
# ------------------------------------------------------------------------------
library(dplyr)       # Data manipulation
library(tidyr)       # Data reshaping (pivot_longer)
library(car)         # Levene's test (leveneTest)
library(ggplot2)     # Residual diagnostic plots
library(openxlsx)    # Excel workbook export

# ------------------------------------------------------------------------------
# 2. Load and prepare data
# ------------------------------------------------------------------------------
data <- read.csv("CWM/RelAbund.csv", sep = ";")
data2 <- data[, 2:20]  # Trait proportion columns (relative abundance)

# Convert relative proportions to percentages (rounded to integers)
data_transformed <- data2 |>
  mutate(across(everything(), ~ round(. * 100, 0)))

data_transformed$BurningRegime <- data$BurningRegime
trait_names <- names(data_transformed)[1:19]

# ------------------------------------------------------------------------------
# 3. Run ANOVA, Tukey, Shapiro-Wilk, and Levene for each trait
# ------------------------------------------------------------------------------

# Storage lists
ovsDB  <- list()   # ANOVA summaries
tuksDB <- list()   # Tukey HSD results
shapDB <- list()   # Shapiro-Wilk normality tests
LevDB  <- list()   # Levene homoscedasticity tests
modDB  <- list()   # ANOVA model objects

for (i in seq_along(trait_names)) {
  x <- trait_names[i]

  # Reshape to long format for a single trait
  dataST <- data_transformed |>
    pivot_longer(cols = all_of(x),
                 names_to  = "Types",
                 values_to = "Total") |>
    select(BurningRegime, Total)

  # One-way ANOVA
  anova_result <- aov(Total ~ BurningRegime, data = dataST)

  # Collect results
  modDB  <- c(modDB,  list(anova_result))
  ovsDB  <- c(ovsDB,  list(summary(anova_result)))
  tuksDB <- c(tuksDB, list(TukeyHSD(anova_result)))
  shapDB <- c(shapDB, list(shapiro.test(residuals(anova_result))))
  LevDB  <- c(LevDB,  list(leveneTest(Total ~ BurningRegime, data = dataST)))
}

# ------------------------------------------------------------------------------
# 4. Print results and diagnostic plots to console
# ------------------------------------------------------------------------------
for (j in seq_along(trait_names)) {
  cat(rep("=", 70), sep = "")
  cat("\n", trait_names[j], "\n")
  cat(rep("=", 70), sep = "")

  cat("\n\nANOVA:\n")
  print(ovsDB[[j]])

  cat("\nTukey HSD post-hoc test:\n")
  print(tuksDB[[j]])

  cat("\nShapiro-Wilk normality test (residuals):\n")
  print(shapDB[[j]])

  cat("\nLevene homogeneity of variance test:\n")
  print(LevDB[[j]])

  # Residuals vs. fitted plot (optional — comment out if not needed)
  modDB2 <- data.frame(
    fitted     = fitted(modDB[[j]]),
    residuals  = residuals(modDB[[j]])
  )
  print(
    ggplot(modDB2, aes(x = fitted, y = residuals)) +
      geom_point() +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      ggtitle(paste("Residuals vs. fitted —", trait_names[j]))
  )

  cat("\n\n")
}

# ------------------------------------------------------------------------------
# 5. Export all results to Excel (one sheet per analysis type)
# ------------------------------------------------------------------------------

# Helper: flatten ANOVA summary to a data.frame row
extract_anova_fg <- function(aov_sum, var_name) {
  df <- as.data.frame(aov_sum[[1]])
  df$Variable <- var_name
  df$Term     <- trimws(rownames(df))
  rownames(df) <- NULL
  df[, c("Variable", "Term", setdiff(names(df), c("Variable", "Term")))]
}

# Helper: flatten TukeyHSD to data.frame
extract_tukey_fg <- function(tuk_obj, var_name) {
  df <- as.data.frame(tuk_obj$BurningRegime)
  df$Variable   <- var_name
  df$Comparison <- rownames(df)
  rownames(df)  <- NULL
  df[, c("Variable", "Comparison",
         setdiff(names(df), c("Variable", "Comparison")))]
}

# Helper: Shapiro-Wilk result
extract_shapiro <- function(shap_obj, var_name) {
  data.frame(
    Variable  = var_name,
    Statistic = shap_obj$statistic,
    p_value   = shap_obj$p.value,
    Normal    = ifelse(shap_obj$p.value > 0.05, "Yes", "No"),
    stringsAsFactors = FALSE
  )
}

# Helper: Levene test result
extract_levene <- function(lev_obj, var_name) {
  df <- as.data.frame(lev_obj)
  data.frame(
    Variable    = var_name,
    F_value     = df[1, "F value"],
    df1         = df[1, "Df"],
    df2         = df[2, "Df"],
    p_value     = df[1, "Pr(>F)"],
    Homogeneous = ifelse(df[1, "Pr(>F)"] > 0.05, "Yes", "No"),
    stringsAsFactors = FALSE
  )
}

# Build each sheet
sheet_anova  <- bind_rows(Map(extract_anova_fg, ovsDB,  trait_names))
sheet_tukey  <- bind_rows(Map(extract_tukey_fg, tuksDB, trait_names))
sheet_shapiro <- bind_rows(Map(extract_shapiro,  shapDB, trait_names))
sheet_levene  <- bind_rows(Map(extract_levene,   LevDB,  trait_names))

# Write workbook
wb <- createWorkbook()

addWorksheet(wb, "ANOVA")
writeData(wb, "ANOVA", sheet_anova)

addWorksheet(wb, "Tukey")
writeData(wb, "Tukey", sheet_tukey)

addWorksheet(wb, "Shapiro")
writeData(wb, "Shapiro", sheet_shapiro)

addWorksheet(wb, "Levene")
writeData(wb, "Levene", sheet_levene)

saveWorkbook(wb, "CWM/ResultsFG.xlsx", overwrite = TRUE)
message("Results exported to CWM/ResultsFG.xlsx")

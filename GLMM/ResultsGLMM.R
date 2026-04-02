# ==============================================================================
# ResultsGLMM.R
# Statistical results — Generalized Linear Mixed Models (GLMM)
# Study: Dung beetle responses (abundance, diversity, removal, functional groups)
#        across burning regimes
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load required libraries
# ------------------------------------------------------------------------------
library(multcomp)    # Simultaneous inference — Tukey contrasts   v1.4-25
library(lme4)        # Mixed-effects models (lmer / glmer.nb)     v1.1-36
library(lmtest)      # Likelihood-ratio tests                      v0.9-40
library(lattice)     # Trellis graphics (diagnostic plots)         v0.21-9
library(tidyverse)   # Data manipulation and visualization         v2.0.0
library(emmeans)     # Estimated marginal means                    v1.10.1
library(psych)       # Descriptive statistics                      v2.4.3
library(performance) # Model diagnostics (overdispersion, residuals) v0.13.0
library(MASS)        # Negative binomial GLMs                      v7.3-60
library(openxlsx)    # Excel workbook export                       (any recent)

# ------------------------------------------------------------------------------
# 2. Load and prepare data
# ------------------------------------------------------------------------------
GLMMDB <- read.csv("GLMM/Data.csv", sep = ";")

GLMMDB$BurningRegime <- factor(GLMMDB$BurningRegime,
                               levels = c("Null", "Low", "High"))

GLMMDB$GrazingEnviroments <- factor(GLMMDB$GrazingEnviroments,
                                    levels = c("SecondaryForest",
                                               "BiodiversePastures",
                                               "GrassMonoculture"))

# ------------------------------------------------------------------------------
# 3. Fit models and collect results
# ------------------------------------------------------------------------------

# Storage lists
aovs    <- list()   # Likelihood-ratio model comparisons (anova null vs full)
models  <- list()   # Model summaries (fixed effects)
ovsDB   <- list()   # Overdispersion check results
resDB   <- list()   # Residual check results
tuksDB  <- list()   # Tukey post-hoc results
var_names <- names(GLMMDB)[4:ncol(GLMMDB)]

# --- 3a. Linear Mixed Models (LMM) for columns 4–7 --------------------------
for (i in 4:7) {
  x <- names(GLMMDB)[i]

  modNull <- as.formula(paste(x, "~ 1 + (1 | GrazingEnviroments)"))
  mod     <- as.formula(paste(x, "~ BurningRegime + (1 | GrazingEnviroments)"))

  DBNull <- lmer(formula = modNull, data = GLMMDB)
  DB     <- lmer(formula = mod,     data = GLMMDB)

  aovs   <- c(aovs,   list(anova(DBNull, DB)))
  tuksDB <- c(tuksDB, list(summary(glht(DB, linfct = mcp(BurningRegime = "Tukey")))))
  models <- c(models, list(summary(DB)))
  ovsDB  <- c(ovsDB,  list(check_overdispersion(DB)))
  resDB  <- c(resDB,  list(check_residuals(DB)))
}

# --- 3b. Negative Binomial GLMMs (GLMM-NB) for columns 8 onward ------------
for (i in 8:ncol(GLMMDB)) {
  x <- names(GLMMDB)[i]

  modNull <- as.formula(paste(x, "~ 1 + (1 | GrazingEnviroments)"))
  mod     <- as.formula(paste(x, "~ BurningRegime + (1 | GrazingEnviroments)"))

  DBNull <- glmer.nb(formula = modNull, data = GLMMDB)
  DB     <- glmer.nb(formula = mod,     data = GLMMDB)

  aovs   <- c(aovs,   list(anova(DBNull, DB)))
  tuksDB <- c(tuksDB, list(summary(glht(DB, linfct = mcp(BurningRegime = "Tukey")))))
  models <- c(models, list(summary(DB)))
  ovsDB  <- c(ovsDB,  list(check_overdispersion(DB)))
  resDB  <- c(resDB,  list(check_residuals(DB)))
}

# ------------------------------------------------------------------------------
# 4. Print results to console
# ------------------------------------------------------------------------------
k <- ncol(GLMMDB) - 3   # Total number of response variables

for (j in seq_len(k)) {
  cat(rep("=", 70), sep = "")
  cat("\n", var_names[j], "\n")
  cat(rep("=", 70), sep = "")
  cat("\n\nModel comparison (LRT):\n")
  print(aovs[[j]])
  cat("\nModel summary:\n")
  print(models[[j]])
  cat("\nTukey post-hoc test:\n")
  print(tuksDB[[j]])
  cat("\nOverdispersion check:\n")
  print(ovsDB[[j]])
  cat("\nResidual check:\n")
  print(resDB[[j]])
  cat("\n\n")
}

# ------------------------------------------------------------------------------
# 5. Export results to Excel (one sheet per analysis type)
# ------------------------------------------------------------------------------

# Helper: flatten a model comparison (anova output) to a data.frame
extract_anova <- function(aov_obj, var_name) {
  df <- as.data.frame(aov_obj)
  df$Variable <- var_name
  df$Model    <- rownames(df)
  rownames(df) <- NULL
  df[, c("Variable", "Model", setdiff(names(df), c("Variable", "Model")))]
}

# Helper: extract fixed-effects table from model summary
extract_summary <- function(sum_obj, var_name) {
  df <- as.data.frame(sum_obj$coefficients)
  df$Variable  <- var_name
  df$Term      <- rownames(df)
  rownames(df) <- NULL
  df[, c("Variable", "Term", setdiff(names(df), c("Variable", "Term")))]
}

# Helper: extract Tukey comparisons
extract_tukey <- function(tuk_obj, var_name) {
  df <- as.data.frame(tuk_obj$test[c("coefficients", "sigma", "tstat", "pvalues")])
  df$Variable   <- var_name
  df$Comparison <- rownames(df)
  rownames(df)  <- NULL
  df[, c("Variable", "Comparison",
         setdiff(names(df), c("Variable", "Comparison")))]
}

# Helper: extract overdispersion result
extract_overdispersion <- function(ov_obj, var_name) {
  data.frame(
    Variable    = var_name,
    Ratio       = if (!is.null(ov_obj$dispersion_ratio))  ov_obj$dispersion_ratio  else NA,
    Chisq       = if (!is.null(ov_obj$chisq_statistic))   ov_obj$chisq_statistic   else NA,
    p_value     = if (!is.null(ov_obj$p_value))           ov_obj$p_value           else NA,
    Overdispersed = if (!is.null(ov_obj$p_value))
      ifelse(ov_obj$p_value < 0.05, "Yes", "No") else NA,
    stringsAsFactors = FALSE
  )
}

# Helper: extract residual check result
extract_residuals <- function(res_obj, var_name) {
  data.frame(
    Variable       = var_name,
    Test           = if (!is.null(res_obj$method))    res_obj$method    else NA,
    Statistic      = if (!is.null(res_obj$statistic)) res_obj$statistic else NA,
    p_value        = if (!is.null(res_obj$p.value))   res_obj$p.value   else NA,
    Interpretation = if (!is.null(res_obj$p.value))
      ifelse(res_obj$p.value > 0.05, "OK (normal)", "Deviation") else NA,
    stringsAsFactors = FALSE
  )
}

# Build each sheet
sheet_anova <- bind_rows(
  Map(extract_anova,         aovs,   var_names)
)
sheet_summary <- bind_rows(
  Map(extract_summary,       models, var_names)
)
sheet_tukey <- bind_rows(
  Map(extract_tukey,         tuksDB, var_names)
)
sheet_overdispersion <- bind_rows(
  Map(extract_overdispersion, ovsDB, var_names)
)
sheet_residuals <- bind_rows(
  Map(extract_residuals,     resDB,  var_names)
)

# Write workbook
wb <- createWorkbook()

addWorksheet(wb, "ANOVA")
writeData(wb, "ANOVA", sheet_anova)

addWorksheet(wb, "Model_Summary")
writeData(wb, "Model_Summary", sheet_summary)

addWorksheet(wb, "Tukey")
writeData(wb, "Tukey", sheet_tukey)

addWorksheet(wb, "Overdispersion")
writeData(wb, "Overdispersion", sheet_overdispersion)

addWorksheet(wb, "Residuals")
writeData(wb, "Residuals", sheet_residuals)

saveWorkbook(wb, "GLMM/ResultsGLMM.xlsx", overwrite = TRUE)
message("Results exported to GLMM/resultsGLMM.xlsx")

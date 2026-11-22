
# ==============================================================================
# SELF-REPORT RATINGS: LINEAR MIXED-EFFECTS ANALYSIS
# ==============================================================================
#
# Description:
#   Compare self-report ratings between external and internal attention
#   conditions using linear mixed-effects models.
#
# Analyses:
#   1. Ease or difficulty of switching attention
#   2. Perceived speed of switching attention
#   3. Attentiveness during attention states
#
# Statistical model:
#   SRR ~ Condition + Site + (1 | Subject)
#
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. INITIALIZATION
# ------------------------------------------------------------------------------

rm(list = ls())

library(lmerTest)
library(emmeans)
library(readxl)
library(car)
library(broom.mixed)


# ------------------------------------------------------------------------------
# 2. INPUT FILE AND ANALYSIS CONFIGURATION
# ------------------------------------------------------------------------------

file_path <- "Scripts/R software for statistical analysis/self_report_ratings.xlsx"

# Names corresponding to the three Excel sheets
analysis_names <- c(
  "Ease or difficulty of switching attention",
  "Perceived speed of switching attention",
  "Attentiveness during attention states"
)

# Condition labels for each analysis
condition_labels <- list(
  c("External switch", "Internal switch"),
  c("External switch", "Internal switch"),
  c("External attention", "Internal attention")
)

# Store results from each analysis
results <- list()


# ==============================================================================
# 3. LINEAR MIXED-EFFECTS ANALYSIS
# ==============================================================================

for (sh in 1:3) {
  
  cat("\n")
  cat("============================================================\n")
  cat("Analysis:", analysis_names[sh], "\n")
  cat("Excel sheet:", sh, "\n")
  cat("============================================================\n")
  
  
  # --------------------------------------------------------------------------
  # 3.1 Load data
  # --------------------------------------------------------------------------
  
  lme_data <- read_excel(
    file_path,
    sheet = sh
  )
  
  str(lme_data)
  
  
  # --------------------------------------------------------------------------
  # 3.2 Prepare categorical variables
  # --------------------------------------------------------------------------
  
  lme_data$Subject <- as.factor(lme_data$Subject)
  lme_data$Condition <- as.factor(lme_data$Condition)
  lme_data$Site <- as.factor(lme_data$Site)
  
  lme_data$Condition <- factor(
    lme_data$Condition,
    levels = c(1, 2),
    labels = condition_labels[[sh]]
  )
  
  
  # --------------------------------------------------------------------------
  # 3.3 Fit linear mixed-effects model
  # --------------------------------------------------------------------------
  
  model1_SRR <- lmer(
    SRR ~ Condition + Site + (1 | Subject),
    data = lme_data
  )
  
  
  # --------------------------------------------------------------------------
  # 3.4 Type III ANOVA
  # --------------------------------------------------------------------------
  
  anova_mdl1_SRR <- anova(
    model1_SRR,
    type = "III"
  )
  
  cat("\nType III ANOVA results:\n")
  print(anova_mdl1_SRR)
  
  
  # --------------------------------------------------------------------------
  # 3.5 Model summary
  # --------------------------------------------------------------------------
  
  cat("\nLinear mixed-effects model summary:\n")
  print(summary(model1_SRR))
  
  
  # --------------------------------------------------------------------------
  # 3.6 Descriptive statistics by condition
  # --------------------------------------------------------------------------
  
  lme_data_ext <- lme_data[
    lme_data$Condition == condition_labels[[sh]][1],
  ]
  
  lme_data_int <- lme_data[
    lme_data$Condition == condition_labels[[sh]][2],
  ]
  
  # External condition
  mean_ext <- mean(lme_data_ext[[3]], na.rm = TRUE)
  sd_ext <- sd(lme_data_ext[[3]], na.rm = TRUE)
  
  # Internal condition
  mean_int <- mean(lme_data_int[[3]], na.rm = TRUE)
  sd_int <- sd(lme_data_int[[3]], na.rm = TRUE)
  
  
  # --------------------------------------------------------------------------
  # 3.7 Display descriptive statistics
  # --------------------------------------------------------------------------
  
  descriptive_stats <- data.frame(
    Condition = condition_labels[[sh]],
    Mean = c(mean_ext, mean_int),
    SD = c(sd_ext, sd_int)
  )
  
  cat("\nDescriptive statistics:\n")
  print(descriptive_stats)
  
  
  # --------------------------------------------------------------------------
  # 3.8 Store results
  # --------------------------------------------------------------------------
  
  results[[sh]] <- list(
    Analysis = analysis_names[sh],
    Model = model1_SRR,
    ANOVA = anova_mdl1_SRR,
    Descriptive_Statistics = descriptive_stats
  )
  
  names(results)[sh] <- paste0("Sheet_", sh)
  
}

# ==============================================================================
# END OF ANALYSIS
# ==============================================================================

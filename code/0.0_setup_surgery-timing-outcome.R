## ================================
## 0. Load needed packages and data
## ================================
# Helper function: install if not already available
load_or_install <- function(pkg){
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# List of needed packages
packages <- c("ggplot2",  "kableExtra", "webshot2","magick","knitr","tidyverse",
              "tidyr","ggh4x","GGally","visdat","car","fastDummies","naniar",
              "rms","Hmisc","patchwork","forcats","xfun","gtsummary","gt","dplyr")
# Load or install each
lapply(packages, load_or_install)
setwd("~/Library/CloudStorage/Dropbox/School/CU/fall 2025/BIOS 6618 adv biostatistical method/final project/surgery_timing_outcome/")
data <- read.csv("./data/Surgery_Timing.csv")

############################################################
## 1. Define variable groups (by role & type)
############################################################

# Primary exposure & outcomes
vars_exposure   <- "hour"
vars_outcomes   <- c("mort30", "complication")

# Demographic
vars_demographic <- c("age", "gender", "race", "bmi")

# Comorbidity indicators (0/1)
vars_comorbidity <- c(
  "baseline_cancer",
  "baseline_cvd",
  "baseline_dementia",
  "baseline_diabetes",
  "baseline_digestive",
  "baseline_osteoart",
  "baseline_psych",
  "baseline_pulmonary"
)

# Risk scores / continuous clinical scores
vars_risk_scores <- c(
  "baseline_charlson",
  "mortality_rsi",
  "complication_rsi",
  "ccsMort30Rate",
  "ccsComplicationRate"
)

# Time-related variables
vars_time <- c("hour", "dow", "month", "moonphase")

# Procedure variable
vars_procedure <- "ahrq_ccs"

# Continuous numeric variables for EDA
vars_continuous <- c(
  "age",
  "bmi",
  "baseline_charlson",
  "mortality_rsi",
  "complication_rsi",
  "ccsMort30Rate",
  "ccsComplicationRate",
  "hour"
)

# Binary variables (0/1)
vars_binary <- c(
  vars_comorbidity,
  vars_outcomes
)

# Categorical variables
vars_categorical <- c(
  "gender",
  "race",
  "asa_status",
  "ahrq_ccs",
  "dow",
  "month",
  "moonphase"
)
n_cols <- ncol(data)
print(length(vars_categorical) + length(vars_binary) + length(vars_continuous) == n_cols)  # should equal ncol(data)
vars_categorical_binary <- c(vars_categorical, vars_binary)
############################################################
## 2. Type conversion
############################################################

data <- data %>%
  # Convert binary variables to factors
  mutate(across(all_of(vars_binary), ~ as.numeric(.))) %>%
  # Convert categorical variables to factors
  mutate(across(all_of(vars_categorical), ~ as.factor(.))) %>%
  # convert continuous variables to numeric
  mutate(across(all_of(vars_continuous), ~ as.numeric(.))
  )

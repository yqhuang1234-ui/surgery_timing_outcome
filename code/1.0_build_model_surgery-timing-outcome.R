#################################
#load data
#################################
# set up the environment
setwd("~/Dropbox/School/CU/fall 2025/BIOS 6618 adv biostatistical method/final project/surgery_timing_outcome")
source("./code/0.0_setup_surgery-timing-outcome.R")
# load processed data
file_name_asa_factor <- "./data/processed/2025-11-24_dropna_regroup-procedure_no-encode.rds"
data <- readRDS(file_name_asa_factor)
str(data)
# check missingness for data
miss_summary <- naniar::miss_var_summary(data)
print(miss_summary)
names(data)

#################################
# build models
#################################
vars_spline <- c("age", "mortality_rsi")
var_exporsure <- "hour"
var_outcome <- "mort30"

# model with 
model_hour_linear <- glm(
  mort30 ~
    rcs(age,4) +                         # 
    gender +
    asa_status +
    baseline_cancer +
    baseline_cvd +
    baseline_psych +
    baseline_pulmonary +
    baseline_charlson +
    rcs(mortality_rsi,4) +               #
    ccsMort30Rate +
    hour +                        # linear for now
    procedure,
  family = binomial,
  data   = data
)

# model with rcs(age, 4)
model_hour_spline <- glm(
  mort30 ~
    rcs(age,4) +                         # 
    gender +
    asa_status +
    baseline_cancer +
    baseline_cvd +
    baseline_psych +
    baseline_pulmonary +
    baseline_charlson +
    rcs(mortality_rsi,4) +               #
    ccsMort30Rate +
    rcs(hour,4) +                        # linear for now
    procedure,
  family = binomial,
  data   = data
)

# lrt _aic/bic comparison
# Likelihood Ratio Test (nested: linear age vs spline age)
anova(model_hour_linear, model_hour_spline, test = "LRT")

# Information criteria
AIC(model_hour_linear, model_hour_spline)
BIC(model_hour_linear, model_hour_spline)

model_core <- glm(
  mort30 ~
    rcs(age, 4) +
    gender +
    asa_status +
    baseline_charlson +
    rcs(mortality_rsi, 4) +
    ccsMort30Rate +
    rcs(hour, 4) +
    procedure,
  family = binomial,
  data   = data
)

model_with_cancer <- update(
  model_core,
  . ~ . + baseline_cancer
)
anova(model_core, model_with_cancer, test = "LRT")

model_with_cvd <- update(
  model_core,
  . ~ . + baseline_cvd
)

anova(model_core, model_with_cvd, test = "LRT")

model_with_psych <- update(
  model_core,
  . ~ . + baseline_psych
)

anova(model_core, model_with_psych, test = "LRT")

model_with_pulm <- update(
  model_core,
  . ~ . + baseline_pulmonary
)

anova(model_core, model_with_pulm, test = "LRT")

summary(model_hour_spline)

anova(model_hour_spline, test="LRT")

library(rms)
plot(Predict(model_hour_spline, hour ))

# 1) Create a sequence of hours
hour_seq <- seq(
  from = min(data$hour, na.rm = TRUE),
  to   = max(data$hour, na.rm = TRUE),
  length.out = 100
)

# 2) Set other covariates to typical values (adjust as needed)
newdat <- data.frame(
  hour              = hour_seq,
  age               = median(data$age, na.rm = TRUE),
  gender            = factor("0", levels = levels(data$gender)),          # pick reference level
  asa_status        = factor("1", levels = levels(data$asa_status)),      # reference ASA
  baseline_cancer   = 0,
  baseline_cvd      = 0,
  baseline_psych    = 0,
  baseline_pulmonary= 0,
  baseline_charlson = median(data$baseline_charlson, na.rm = TRUE),
  mortality_rsi     = median(data$mortality_rsi, na.rm = TRUE),
  ccsMort30Rate     = median(data$ccsMort30Rate, na.rm = TRUE),
  procedure         = factor("Orthopedic & Spine", levels = levels(data$procedure))  # choose ref
)

# 3) Get predicted log-odds + SE
pred <- predict(model_hour_spline, newdata = newdat, type = "link", se.fit = TRUE)

# 4) Convert to probabilities + CI
newdat$pred_prob  <- plogis(pred$fit)
newdat$lower_ci   <- plogis(pred$fit - 1.96 * pred$se.fit)
newdat$upper_ci   <- plogis(pred$fit + 1.96 * pred$se.fit)

# 5) Plot
ggplot(newdat, aes(x = hour, y = pred_prob)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2) +
  labs(
    x = "Surgery start time (hour)",
    y = "Predicted 30-day mortality",
    title = "Adjusted relationship between surgery time and 30-day mortality"
  ) +
  theme_minimal()

#fit model with categorical 
# Create binary hour category
data$hour_bin <- ifelse(data$hour < 12, "early", "late")
data$hour_bin <- factor(data$hour_bin, levels = c("early", "late"))
model_hour_bin <- glm(
  mort30 ~ 
    rcs(age, 4) +
    gender +
    asa_status +
    baseline_cancer +
    baseline_cvd +
    baseline_psych +
    baseline_pulmonary +
    baseline_charlson +
    rcs(mortality_rsi, 4) +
    ccsMort30Rate +
    hour_bin +
    procedure,
  data   = data,
  family = binomial
)

AIC(model_hour_linear, model_hour_spline, model_hour_bin)
BIC(model_hour_linear, model_hour_spline, model_hour_bin)
summary(model_hour_bin)
anova(model_hour_bin, test="LRT")

data$hour_quartile <- cut(
  data$hour,
  breaks = quantile(data$hour, probs = c(0, 0.25, 0.5, 0.75, 1),
                    na.rm = TRUE),
  include.lowest = TRUE,
  labels = c("Q1", "Q2", "Q3", "Q4")
)
model_hour_quartile <- glm(
  mort30 ~ 
    rcs(age, 4) +
    gender +
    asa_status +
    baseline_cancer +
    baseline_cvd +
    baseline_psych +
    baseline_pulmonary +
    baseline_charlson +
    rcs(mortality_rsi, 4) +
    ccsMort30Rate +
    hour_quartile +
    procedure,
  data   = data,
  family = binomial
)
AIC(model_hour_linear, model_hour_spline, model_hour_bin, model_hour_quartile)
BIC(model_hour_linear, model_hour_spline, model_hour_bin, model_hour_quartile)
anova(model_hour_quartile, test="LRT")


data$hour_group <- cut(
  data$hour,
  breaks = c(7, 10, 13, 16, 19),
  include.lowest = TRUE,
  labels = c("7–10", "10–13", "13–16", "16–19")
)
model_hour_group <- glm(
  mort30 ~ 
    rcs(age, 4) +
    gender +
    asa_status +
    baseline_cancer +
    baseline_cvd +
    baseline_psych +
    baseline_pulmonary +
    baseline_charlson +
    rcs(mortality_rsi, 4) +
    ccsMort30Rate +
    hour_group +
    procedure,
  data   = data,
  family = binomial
)
AIC(model_hour_linear, model_hour_spline, model_hour_bin, model_hour_quartile, model_hour_group)
BIC(model_hour_linear, model_hour_spline, model_hour_bin, model_hour_quartile, model_hour_group)
anova(model_hour_group, test="LRT")

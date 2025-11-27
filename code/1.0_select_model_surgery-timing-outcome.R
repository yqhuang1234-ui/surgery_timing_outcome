#################################
#load data
#################################
# set up the environment
setwd("~/Dropbox/School/CU/fall 2025/BIOS 6618 adv biostatistical method/final project/surgery_timing_outcome")
source("./code/0.0_setup_surgery-timing-outcome.R")
# load processed data
file_name_asa_factor <- "./data/processed/2025-11-25_dropna_regroup-procedure_no-encode.rds"
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
#get aic and bic
AIC(model_hour_spline)
BIC(model_hour_spline)

model_hour_spline_linear <- glm(
  mort30 ~
    age +                         # 
    gender +
    asa_status +
    baseline_cancer +
    baseline_cvd +
    baseline_psych +
    baseline_pulmonary +
    baseline_charlson +
    rcs(mortality_rsi,3) +               #
    ccsMort30Rate +
    hour +                        # linear for now
    procedure,
  family = binomial,
  data   = data
)
summary(model_hour_spline_linear)
AIC(model_hour_spline_linear)
BIC(model_hour_spline_linear)
summary(model_hour_spline)
anova(model_hour_spline, test="LRT")
model_hour_spline_complication <- glm(
  complication ~
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

summary(model_hour_spline_complication)
anova(model_hour_spline_complication, test="LRT")

# count by hour

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

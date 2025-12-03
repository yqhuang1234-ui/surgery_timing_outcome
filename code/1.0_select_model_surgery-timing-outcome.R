#################################
#load data
#################################
# set up the environment
setwd("~/Dropbox/School/CU/fall 2025/BIOS 6618 adv biostatistical method/final project/surgery_timing_outcome")
source("./code/0.0_setup_surgery-timing-outcome.R")
# load processed data
file_name_asa_factor <- "./data/processed/2025-12-01_dropna_regroup-procedure_no-encode.rds"
data <- readRDS(file_name_asa_factor)
data$hour_cat <- relevel(data$hour_cat, ref = "8-11am")


str(data)
# check missingness for data
miss_summary <- naniar::miss_var_summary(data)
print(miss_summary)
names(data)
dim(data)


#################################
# build models
#################################
# set up
vars_spline <- c("age", "mortality_rsi")
var_exposure <- "hour"
var_exposure_cat <- "hour_cat"
var_outcome <- "mort30"
dd <- datadist(data)
options(datadist = "dd")

# automatically define other covariates
other_covariates <- setdiff(
  names(data),
  c(var_outcome, var_exposure, vars_spline, var_exposure_cat)
)
print(other_covariates)
print(length(other_covariates))

# construct base rhs
base_rhs <- paste(
  "rcs(age, K) + rcs(mortality_rsi, K) + rcs(hour, K)",
  if (length(other_covariates) > 0)
    paste("+", paste(other_covariates, collapse = " + "))
  else ""
)

# fit models with different K
fit_k <- list()
AIC_k <- numeric()

for (K in c(0, 3, 4,5)) {
  if (K == 0) {
    # k = 0 = linear model: no rcs(), just linear terms
    f <- lrm(
      as.formula(
        paste("mort30 ~ age + mortality_rsi + hour",
              if (length(other_covariates) > 0)
                paste("+", paste(other_covariates, collapse = " + "))
              else "")
      ),
      data = data
    )
  } else {
    f <- lrm(
      as.formula(paste("mort30 ~", gsub("K", K, base_rhs))),
      data = data
    )
  }
  fit_k[[as.character(K)]] <- f
  AIC_k[as.character(K)] <- round(AIC(f),0)
}

AIC_k
# choose the k with lowest AIC
best_k <- as.numeric(names(which.min(AIC_k)))
print(paste("Best K:", best_k))
best_model <- fit_k[[as.character(best_k)]]

# within best k, check which variables are truely nonlinear
#if nonlinear has small p keep spline, if p large change to linear
anova(best_model)

# -------------------------
# FIT MODEL WITH HOUR AS CATEGORICAL VARIABLE
# -------------------------
final_model_rhs <- as.formula(
  paste("mort30 ~ age + rcs(mortality_rsi,3) + hour_cat",
        if (length(other_covariates) > 0)
          paste("+", paste(other_covariates, collapse = " + "))
        else "")
)
fit_final_basic <- lrm(
  final_model_rhs,
  data = data
)
anova(fit_final_basic)

#check AIC for all models
c(
  AIC_k,
  final=round(AIC(fit_final_basic),0)
)

# unadjusted effect of hour
fit_unadjusted <- lrm(
  mort30 ~ hour_cat,
  data = data
)
plot(Predict(fit_unadjusted, hour_cat, fun=plogis), xlab="Surgery Hour", ylab="Predicted 30-day Mortality Probability", main="Unadjusted Effect of Surgery Hour on 30-day Mortality")
plot(Predict(fit_final_basic, hour_cat))

#################################
# check overly influential points
#################################
lm_fit2 <- glm(final_model_rhs, data=data, family=binomial)
summary(lm_fit2)
# Jackknife
jackknife_outliers <- sum(abs(rstudent(lm_fit2)) > 3)

# Cook's D
cooks_d_outliers <- sum(cooks.distance(lm_fit2) > 1)

# Leverage
leverage_threshold <- 2 * (length(coef(lm_fit2)) / nrow(data))
leverage_outliers <- sum(hatvalues(lm_fit2) > leverage_threshold)

# DFBETAS Intercept
dfbetas_threshold <- 2 / sqrt(nrow(data))
dfbetas_intercept_outliers <- sum(abs(dfbetas(lm_fit2)[,1]) > dfbetas_threshold)
# DFBETAS Slope
dfbetas_slope_outliers <- sum(abs(dfbetas(lm_fit2)[,2]) > dfbetas_threshold)

# DFFITS
dffits_threshold <- 2 * sqrt(length(coef(lm_fit2)) / nrow(data))
dffits_outliers <- sum(abs(dffits(lm_fit2)) > dffits_threshold)

# Create table
outlier_table <- data.frame(
  Method = c("Jackknife ±3", "Cook's D", "Leverage", "DFBETAS Intercept", "DFBETAS Slope", "DFFITS"),
  Number_of_Observations = c(jackknife_outliers, cooks_d_outliers, leverage_outliers, dfbetas_intercept_outliers, dfbetas_slope_outliers, dffits_outliers)
)
kable(outlier_table)

###################################################
# plots for report purposes
###################################################
#save_plot
# create figure directory if not exist
today_date <- format(Sys.Date(), "%Y-%m-%d")
figure_dir <- file.path("results", today_date,'/')
print(figure_dir)

if (!dir.exists(figure_dir)) {
  dir.create(figure_dir, recursive = TRUE)
}

#################################
# get results from final model
#################################
# this object contains coeffients, se, z statistic, wald p-values, confidence intervals
## 1) Coefficients (betas)
co <- coef(fit_final_basic)

## 2) Variance–covariance matrix and SEs
vc <- vcov(fit_final_basic)
se <- sqrt(diag(vc))

## 3) Wald Z-statistics and p-values
z  <- co / se
p  <- 2 * pnorm(-abs(z))

## 4) Wald 95% CI for betas
ci_beta <- confint.default(fit_final_basic)  # works for lrm objects

# -------------- VARIABLE LABELS & LEVEL LABELS -------------------

## 5) Put everything together on OR scale
options(scipen = 999)   # prevents scientific notation
or_table <- data.frame(
  term    = names(co),
  beta    = round(co, 3),
  se      = round(se, 3),
  z       = round(z, 3),
  p_value = ifelse(p < 0.001, "<0.001", round(p, 3)),
  OR      = round(exp(co), 2),
  lower95 = round(exp(ci_beta[, 1]), 2),
  upper95 = round(exp(ci_beta[, 2]), 2),
  row.names = NULL,
  check.names = FALSE
)
or_table <- or_table %>%
  mutate(
    clean_term = dplyr::recode(term,
                               
                               # continuous variables  
                               "age"               = "Age",
                               "mortality_rsi"     = "Risk Stratification Index (30-day mortality)",
                               "baseline_charlson" = "Charlson Comorbidity Index",
                               
                               # gender (binary)
                               "gender=1" = "Gender=Female",
                               
                               # ASA levels
                               "asa_status=2" = "ASA status III",
                               "asa_status=3" = "ASA status IV–VI",
                               
                               # Procedure categories
                               "procedure=GI"     = "Gastrointestinal & Hernia (vs other)",
                               "procedure=Ortho"  = "Orthopedic & Spine (vs other)",
                               "procedure=Gyn"    = "Gynecologic (vs other)",
                               "procedure=Uro"    = "Urologic (vs other)",
                  
                               
                               # Binary comorbidities (0/1)
                               "baseline_cancer"    = "Cancer",
                               "baseline_cvd"       = "Cardiovascular/Cerebrovascular disease",
                               "baseline_dementia"  = "Dementia",
                               "baseline_diabetes"  = "Diabetes",
                               "baseline_digestive" = "Digestive disease",
                               "baseline_osteoart"  = "Osteoarthritis",
                               "baseline_psych"     = "Psychiatric disorder",
                               "baseline_pulmonary" = "Pulmonary disease",
                               
                               # Hour categories
                               "hour_cat=8-11am"   = "Surgery Hour 8–11am",
                               "hour_cat=11-1pm" = "Surgery Hour 11–1pm",
                               "hour_cat=1-3pm"    = "Surgery Hour 1–3pm",
                               "hour_cat=3-5pm"    = "Surgery Hour 3–5pm",
                               "hour_cat=5-7pm"    = "Surgery Hour 5–7pm",
                               "hour_cat=6-8am"   = "Surgery Hour 6–8am",
                               
                               .default = term    # fallback preserves anything unmapped
    )
  )


#significant variables at alpha = 0.05
or_table %>%
  filter(p_value != "<0.001" & p_value < 0.05 | p_value == "<0.001")

plot_df <- or_table %>%
  # remove vector of terms don't want to plot
  filter(!term %in% c("Intercept","mortality_rsi","mortality_rsi'")) %>%
  filter(!is.na(OR)) %>%
  mutate(
    # numeric p-value
    p_num = 2 * pnorm(-abs(beta / se)),
    
    # significance flag
    signif = p_num < 0.05 &
      ((lower95 > 1) | (upper95 < 1)),
    
    signif_group = ifelse(
      signif,
      "Significant (p < 0.05; 95% CI does not include 1)",
      "Not significant (p ≥ 0.05 or 95% CI includes 1)"
    ),
    
    # formatted p-value
    label_p = ifelse(
      p_num < 0.001, 
      "p<0.001",
      paste0("p=", sprintf("%.3f", p_num))
    ),
    
    # text label: ASA IV–VI gets 2 lines, others stay 1 line
    label = if_else(
      clean_term == "ASA status IV–VI",
      paste0(sprintf("%.2f", OR), "\n(", label_p, ")"),
      paste0(sprintf("%.2f", OR), " (", label_p, ")")
    )
  ) %>%
  arrange(OR) %>%
  mutate(clean_term = factor(clean_term, levels = clean_term))

plot_df <- plot_df %>%
  mutate(
    OR_text = sprintf("%.2f", OR),
    p_text  = paste0("(", label_p, ")"),
    v_OR = if_else(clean_term == "ASA status IV–VI",  0.3, 0),
    v_p  = if_else(clean_term == "ASA status IV–VI", -0.3, 0)
  )


ggplot(plot_df, aes(x = OR, y = clean_term, color = signif_group)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  
  geom_errorbarh(aes(xmin = lower95, xmax = upper95), height = 0.15) +
  geom_point(size = 2) +
  
  geom_text(
    aes(x = upper95, label = label),
    hjust = -0.03,
    size = 3
  ) +
  
  scale_x_log10(expand = expansion(mult = c(0.05, 0.20))) +
  
  scale_color_manual(
    name = NULL,
    values = c(
      "Significant (p < 0.05; 95% CI does not include 1)" = "#CC3311",
      "Not significant (p ≥ 0.05 or 95% CI includes 1)"    = "grey40"
    ),
    labels = c("Not significant", "Significant"),
    breaks = c(
      "Not significant (p ≥ 0.05 or 95% CI includes 1)",
      "Significant (p < 0.05; 95% CI does not include 1)"
    )
  ) +
  
  guides(color = guide_legend(nrow = 1)) +
  
  labs(
    x = "Odds ratio (log scale)",
    y = NULL
  ) +
  
  theme_light() +
  theme(
    panel.background  = element_rect(fill = "white", color = NA),
    panel.border      = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor.y = element_blank(),
    
    axis.line   = element_blank(),
    axis.ticks  = element_blank(),
    # x title no longer bold
    axis.title.x = element_text(size = 12, face = "plain"),
    axis.text   = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, hjust = 0),
    
    legend.position = "top",
    legend.title    = element_blank(),
    legend.text     = element_text(size = 9),
    legend.key      = element_blank(),
    
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin     = margin(2, 2, 2, 2)
  )

#save figure
ggsave(paste0(figure_dir, "2_final-model-OR-figure.png"), width=5.6, height=5.5, dpi = 600)

# effect plot
# plot unadjusted and adjusted effects of hour on mortality showing two lines


# Use fun = plogis to convert log-odds -> probability (0–1)
# Get predictions for unadjusted and adjusted models
p1 <- as.data.frame(Predict(fit_unadjusted, hour_cat, fun = plogis)) %>%
  mutate(model = "Unadjusted")

p2 <- as.data.frame(Predict(fit_final_basic, hour_cat, fun = plogis)) %>%
  mutate(model = "Adjusted")

combined_plot <- bind_rows(p1, p2)
combined_plot$hour_cat_plot <- factor(
  combined_plot$hour_cat,
  levels = c("6-8am", "8-11am", "11-1pm", "1-3pm", "3-5pm", "5-7pm")
)
library(scales)  # for label_percent()

# Ensure hour_cat is a factor in the combined data too

ggplot(
  combined_plot,
  aes(x = hour_cat_plot, y = yhat, color = model, group = model)
) +
  # 95% CI ribbons
  geom_ribbon(
    aes(ymin = lower, ymax = upper, fill = model),
    alpha = 0.20,
    color = NA
  ) +
  # Lines for predicted probabilities
  geom_line(size = 1.3) +
  
  labs(
    x = "Surgery hour category",
    y = "Predicted 30-day mortality (%)"
  ) +
  
  # Colors for unadjusted vs adjusted
  scale_color_manual(
    name = NULL,
    values = c(
      "Unadjusted" = "#0BA6DF",
      "Adjusted"   = "#CC3311"
    )
  ) +
  scale_fill_manual(
    name = NULL,
    values = c(
      "Unadjusted" = "#0BA6DF",
      "Adjusted"   = "#CC3311"
    )
  ) +
  
  # Show y-axis as percentages (but keep underlying data as probabilities)
  scale_y_continuous(
    labels = label_percent(accuracy = 0.01),  # e.g. 0.45%
    expand = expansion(mult = c(0, 0.05))
  ) +
  
  # JAMA-ish theme
  theme_bw() +
  theme(
    panel.background   = element_rect(fill = "white", color = NA),
    panel.border       = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.4),
    
    axis.line    = element_blank(),
    axis.ticks   = element_blank(),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x  = element_text(size = 10),
    axis.text.y  = element_text(size = 10),
    
    legend.position = "top",
    legend.title    = element_blank(),
    legend.key      = element_blank(),
    legend.text     = element_text(size = 10),
    
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin     = margin(5, 5, 5, 5)
  )
#1A4E99
# save figure
ggsave(paste0(figure_dir, "2_hour-effect-compare-figure.png"), width=5, height=4, dpi = 600)

#effect plot mortality_rsi
df_rsi <- as.data.frame(
  Predict(fit_final_basic, mortality_rsi, fun = plogis)
) %>%
  mutate(
    mort_pct = yhat * 100,
    lower_pct = lower * 100,
    upper_pct = upper * 100
  )

ggplot(df_rsi, aes(x = mortality_rsi, y = mort_pct)) +
  geom_line(size = 1.2, color = "#0BA6DF") +
  geom_ribbon(aes(ymin = lower_pct, ymax = upper_pct),
              alpha = 0.20, fill = "#0BA6DF") +
  labs(
    x = "Risk Stratification Index (30-day Mortality)",
    y = "Predicted 30-day Mortality (%)"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

# save figure
ggsave(paste0(figure_dir, "2_mortality-rsi-effect_figure.png"), width=4, height=3, dpi = 600)

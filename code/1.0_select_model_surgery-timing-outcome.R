#################################
#load data
#################################
# set up the environment
setwd("~/Dropbox/School/CU/fall 2025/BIOS 6618 adv biostatistical method/final project/surgery_timing_outcome")
source("./code/0.0_setup_surgery-timing-outcome.R")
# load processed data
file_name_asa_factor <- "./data/processed/2025-11-30_dropna_regroup-procedure_no-encode.rds"
data <- readRDS(file_name_asa_factor)

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
var_exporsure <- "hour"
var_exporsure_cat <- "hour_cat"
var_outcome <- "mort30"
dd <- datadist(data)
options(datadist = "dd")

# automatically define other covariates
other_covariates <- setdiff(
  names(data),
  c(var_outcome, var_exporsure, vars_spline, var_exporsure_cat)
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
# get results
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


or_table

#significant variables at alpha = 0.05
or_table %>%
  filter(p_value != "<0.001" & p_value < 0.05 | p_value == "<0.001")

# 1. Start from your existing or_table
# or_table has: term, OR, lower95, upper95, p_value, etc.



  
plot_df <- or_table %>%
    filter(term != "Intercept") %>%
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
      
      # format p-value
      label_p = ifelse(
        p_num < 0.001, 
        "p<0.001",
        paste0("p=", sprintf("%.3f", p_num))
      ),
      
      # text label for OR
      label = paste0(sprintf("%.2f", OR), " (", label_p, ")")
    ) %>%
    # sort descending by OR
    arrange(OR) %>%
    mutate(term = factor(term, levels = term))

ggplot(plot_df, aes(x = OR, y = term, color = signif_group)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  
  geom_errorbarh(aes(xmin = lower95, xmax = upper95), height = 0.15) +
  geom_point(size = 2) +
  
  # label placed at upper CI bound
  geom_text(
    aes(x = upper95, label = label),
    hjust = -0.1,
    size = 3
  ) +
  
  # log OR axis
  scale_x_log10(expand = expansion(mult = c(0.05, 0.40))) +
  
  scale_color_manual(
    name = NULL,
    values = c(
      "Significant (p < 0.05; 95% CI does not include 1)" = "firebrick",
      "Not significant (p ≥ 0.05 or 95% CI includes 1)" = "grey40"
    )
  ) +
  guides(color = guide_legend(ncol = 2)) + 
  labs(
    x = "Odds ratio (log scale)",
    y = NULL
  ) +
  
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 9),
    legend.position = "top",
    legend.justification = "left"
  )
#save_plot
#add today's date to file name
today_date <- format(Sys.Date(), "%Y-%m-%d")
print(today_date)
#create directory if not exist
if (!dir.exists("./data/processed/")) {
  dir.create("./data/processed/")
}
#save data
filename_without_encode <- paste0(today_date, "_dropna_regroup-procedure_no-encode.rds")
saveRDS(data_to_model, file = paste0("./data/processed/", filename_without_encode))
cat("Cleaned data saved for modeling:\n")


# effect plot
# plot unadjusted and adjusted effects of hour on mortality showing two lines


# Use fun = plogis to convert log-odds -> probability (0–1)
# Get predictions for unadjusted and adjusted models
p1 <- as.data.frame(Predict(fit_unadjusted, hour_cat, fun = plogis)) %>%
  mutate(model = "Unadjusted")

p2 <- as.data.frame(Predict(fit_final_basic, hour_cat, fun = plogis)) %>%
  mutate(model = "Adjusted")

combined <- bind_rows(p1, p2)

# Ensure hour_cat is a factor in the combined data too

ggplot(
  combined,
  aes(x = hour_cat, y = yhat, color = model, group = model)
) +
  geom_line(aes(group = model), size = 1.2) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper, fill = model, group = model),
    alpha = 0.2,
    color = NA
  ) +
  labs(
    x = "Surgery Hour Category",
    y = "Predicted 30-day Mortality Probability"
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.title = element_blank()   # remove legend title
  )
  scale_color_manual(values = c("Unadjusted" = "red", "Adjusted" = "blue")) +
  scale_fill_manual(values = c("Unadjusted" = "red", "Adjusted" = "blue"))

#effect plot mortality_rsi
ggplot(as.data.frame(Predict(fit_final_basic, mortality_rsi, fun=plogis)), 
       aes(x = mortality_rsi, y = yhat)) +
  geom_line(size=1.2, color="blue") +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, fill="blue") +
  labs(
    x = "Mortality RSI",
    y = "Predicted 30-day Mortality Probability",
    title = "Adjusted Association Between Mortality RSI and 30-day Mortality"
  ) +
  theme_bw()

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



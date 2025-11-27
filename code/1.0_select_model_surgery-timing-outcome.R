#################################
#load data
#################################
# set up the environment
setwd("~/Dropbox/School/CU/fall 2025/BIOS 6618 adv biostatistical method/final project/surgery_timing_outcome")
source("./code/0.0_setup_surgery-timing-outcome.R")
# load processed data
file_name_asa_factor <- "./data/processed/2025-11-27_dropna_regroup-procedure_no-encode.rds"
data <- readRDS(file_name_asa_factor)
str(data)
# check missingness for data
miss_summary <- naniar::miss_var_summary(data)
print(miss_summary)
names(data)

#################################
# build models
#################################
# set up
vars_spline <- c("age", "mortality_rsi")
var_exporsure <- "hour"
var_outcome <- "mort30"
dd <- datadist(data)
options(datadist = "dd")

# automatically define other covariates
other_covariates <- setdiff(
  names(data),
  c(var_outcome, var_exporsure, vars_spline)
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

for (K in c(0, 3, 4)) {
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
  AIC_k[as.character(K)] <- AIC(f)
}

AIC_k
# choose the k with lowest AIC
best_k <- as.numeric(names(which.min(AIC_k)))
print(paste("Best K:", best_k))
best_model <- fit_k[[as.character(best_k)]]

# within best k, check which variables are truely nonlinear
#if nonlinear has small p keep spline, if p large change to linear
anova(best_model)

# finalize model formula
fit_final_basic <- lrm(
  as.formula(
    paste("mort30 ~ age + rcs(mortality_rsi,3) + hour",
          if (length(other_covariates) > 0)
            paste("+", paste(other_covariates, collapse = " + "))
          else "")
  ),
  data = data
)
anova(fit_final_basic)

#check AIC for all models
c(
  AIC_k,
  final=AIC(fit_final_basic)
)

plot(Predict(fit_final_basic, hour, fun=plogis), xlab="Surgery Hour", ylab="Predicted 30-day Mortality Probability", main="Effect of Surgery Hour on 30-day Mortality")
plot(Predict(fit_final_basic, mortality_rsi, fun=plogis), xlab="Mortality RSI", ylab="Predicted 30-day Mortality Probability", main="Effect of Mortality RSI on 30-day Mortality")
plot(Predict(fit_final_basic, age, fun=plogis), xlab="Age", ylab="Predicted 30-day Mortality Probability", main="Effect of Age on 30-day Mortality")

# unadjusted effect of hour
fit_unadjusted <- lrm(
  mort30 ~ hour,
  data = data
)
plot(Predict(fit_unadjusted, hour, fun=plogis), xlab="Surgery Hour", ylab="Predicted 30-day Mortality Probability", main="Unadjusted Effect of Surgery Hour on 30-day Mortality")

#################################
# model diagnostics
#################################
stats <- fit_final_basic$stats
ggplot( Predict(fit_final_basic), sepdiscrete= 'vertical' , vnames = 'names' ,
         rdata =data ,
         histSpike.opts = list(frac= function(fit_final_basic) .1*fit_final_basic/max(fit_final_basic) ))



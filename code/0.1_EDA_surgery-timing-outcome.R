## ================================
## 0. Load needed packages and data
## ================================
source("./code/0.0_setup_surgery-timing-outcome.R")
install.packages(c("visdat"))
library(visdat)
head(data)
colnames(data)
data <- data %>%
  # Convert binary variables to factors
  mutate(across(all_of(vars_binary), ~ as.numeric(.))) %>%
  # Convert categorical variables to factors
  mutate(across(all_of(vars_categorical), ~ as.factor(.))) %>%
  # convert continuous variables to numeric
  mutate(across(all_of(vars_continuous), ~ as.numeric(.))
  )
############################################################
## 3. Quick data overview & missingness
############################################################

cat("=== Structure ===\n")
str(data)

cat("\n=== Basic summary ===\n")
summary(data)

cat("\n=== Missingness by variable (proportion) ===\n")
miss_prop <- colMeans(is.na(data))
miss_vars <- names(miss_prop[miss_prop > 0])
print(round(miss_prop[miss_prop >0] * 100, 4))
# summarize missingness by outcome groups
for (outcome in vars_outcomes) {
  cat(paste0("\n=== Missingness by variable (proportion) stratified by ", outcome, " ===\n"))
  miss_by_outcome <- data %>%
    group_by(!!sym(outcome)) %>%
    summarise(across(everything(), ~ mean(is.na(.)))) %>%
    select(all_of(miss_vars))
  print(round(miss_by_outcome * 100, 4))
}

############################################################
## Pairwise plots for selected continuous variables
############################################################
plot_corr <- GGally::ggcorr(
  data %>% select(all_of(vars_continuous)),
  label = TRUE,
  label_alpha = TRUE,
  hjust = 0.75,
  size = 3,
  low = "blue",
  mid = "white",
  high = "red",
  nbreaks = 6
) +
  ggtitle("Correlation matrix of continuous variables") +
  theme_minimal()
print(plot_corr)
# scatter plot mortality_rsi vs ccsMort30Rate 
ggplot(data, aes(x = mortality_rsi, y = ccsMort30Rate)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    title = "Scatter plot of Mortality RSI vs CCS Mort30 Rate",
    x = "Mortality RSI",
    y = "CCS Mort30 Rate"
  ) +
  theme_minimal()
############################################################
## Base R EDA for continuous variables (subplots via par)
############################################################

print(length(vars_continuous))
# Histograms
par(mfrow = c(3, 3))
for (v in vars_continuous) {
  hist(data[[v]],
       main = paste("Histogram of", v),
       xlab = v,
       col = "grey",
       border = "white")
}
# Boxplots to see outliers
par(mfrow = c(3, 3))
for (v in vars_continuous) {
  boxplot(data[[v]],
          main = paste("Boxplot of", v),
          ylab = v,
          col = "lightblue")
}


############################################################
## Base R EDA for categorical variables (subplots)
############################################################

cat_vars <- vars_categorical + vars_binary
print(length(cat_vars))
par(mfrow = c(4, 3))
for (v in cat_vars[1:12]) {
  tab <- table(data[[v]])
  barplot(tab,
          main = paste("Barplot of", v),
          ylab = "Count",
          las = 2)
}
group_by(data, mort30) %>%
  summarise(n = n())
group_by(data, complication) %>%
  summarise(n = n())

par(mfrow = c(3, 2))
for (v in cat_vars[13:length(cat_vars)]) {
  tab <- table(data[[v]])
  barplot(tab,
          main = paste("Barplot of", v),
          ylab = "Count",
          las = 2)
}

############################################################
## Relationships: exposure (hour) vs outcome (mort30)
############################################################

# histogram of hour
ggplot(data, aes(x = hour)) +
  geom_histogram(binwidth = 1, fill = "lightblue", color = "black") +
  labs(
    title = "Histogram of Hour of Surgery",
    x = "Hour of Surgery",
    y = "Count"
  ) +
  theme_minimal()

# Crude mortality by hour (binned)
data %>%
  mutate(hour_bin = cut(hour, breaks = seq(0, 24, by = 2), right = FALSE)) %>%
  group_by(hour_bin) %>%
  summarise(
    mort_rate = mean(mort30, na.rm = TRUE),
    n = n()
  ) %>%
  ggplot(aes(x = hour_bin, y = mort_rate, group = 1)) +
  geom_point() +
  geom_line() +
  labs(
    title = "Crude 30-day mortality rate by hour of surgery (2-hour bins)",
    x = "Hour of surgery (binned)",
    y = "Mortality rate"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Smooth curve (loess) of mort30 vs hour
ggplot(data, aes(x = hour, y = mort30)) +
  geom_smooth(method = "loess", se = TRUE) +
  labs(
    title = "Smoothed relationship between hour of surgery and 30-day mortality",
    x = "Hour of surgery",
    y = "Mortality rate"
  ) +
  theme_minimal()
############################################################
## relationship between mortality rate and categorical variables
############################################################

# Mort30 and asa status
data %>%
  dplyr::group_by(ahrq_ccs) %>%
  dplyr::summarise(
    mort_rate = mean(mort30, na.rm = TRUE),
    n = n()
  ) %>%
  ggplot(aes(x = ahrq_ccs, y = mort_rate)) +
  geom_col() +
  geom_text(
    aes(label = sprintf("%.3f", mort_rate)),
    vjust = -0.3,
    size  = 3
  ) +
  labs(
    title = "30-day mortality proportion by ahrq_ccs",
    x     = "AHRQ CCS",
    y     = "Mortality proportion"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1) 
  )
unique(data$mort30)
# show # of mortality cases by ahrq_ccs
data %>%
  dplyr::group_by(ahrq_ccs) %>%
  dplyr::summarise(
    mort_cases = sum(mort30, na.rm = TRUE),
    n          = n(),
    .groups    = "drop"
  ) 
# histogram of ahrq_ccs and show mortality rate as text
data %>%
  dplyr::group_by(ahrq_ccs) %>%
  dplyr::summarise(
    mort_rate = mean(mort30, na.rm = TRUE),
    n         = n(),
    .groups   = "drop"
  ) %>%
  ggplot(aes(x = ahrq_ccs, y = n)) +
  geom_col() +   # use the precomputed counts
  geom_text(
    aes(
      y     = n + 0.03 * max(n),             # a bit above the bar
      label = sprintf("%.3f", mort_rate)
    ),
    size = 3
  ) +
  labs(
    title = "Counts of AHRQ CCS with 30-day Mortality Rate",
    x     = "AHRQ CCS",
    y     = "Count"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


# mortality_rsi by asa status
ggplot(data, aes(x = asa_status, y = mortality_rsi)) +
  geom_boxplot() +
  labs(
    title = "Boxplot of Mortality RSI by ASA Status",
    x = "ASA Status",
    y = "Mortality RSI"
  ) +
  theme_minimal()

plot_mortality_bar <- function(data, cat_var) {
  
  var_label <- rlang::as_label(rlang::enquo(cat_var))
  
  data %>%
    group_by({{ cat_var }}) %>%
    summarise(
      mort_rate = mean(mort30, na.rm = TRUE),
      n         = n(),
      .groups   = "drop"
    ) %>%
    ggplot(aes(x = {{ cat_var }}, y = mort_rate)) +
    geom_col() +
    geom_text(
      aes(label = sprintf("%.3f", mort_rate)),
      vjust = -0.3,
      size  = 3
    ) +
    labs(
      title = paste("30-day mortality proportion by", var_label),
      x     = var_label,
      y     = "Mortality proportion"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}
plot_mortality_bar(data, complication)
plot_mortality_bar(data, asa_status)
plot_rsi_box <- function(data, cat_var) {
  
  var_label <- rlang::as_label(rlang::enquo(cat_var))
  
  ggplot(data, aes(x = {{ cat_var }}, y = mortality_rsi)) +
    geom_boxplot() +
    labs(
      title = paste("Boxplot of Mortality RSI by", var_label),
      x     = var_label,
      y     = "Mortality RSI"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}
plot_rsi_box(data, asa_status)


# mortality_rsi by ahrq_ccs
mort_rates <- data %>%
  group_by(ahrq_ccs) %>%
  summarise(
    mort_rate = mean(mort30, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(data, aes(x = ahrq_ccs, y = mortality_rsi)) +
  geom_boxplot() +
  geom_text(
    data = mort_rates,
    aes(
      x = ahrq_ccs,
      y = Inf,
      label = sprintf("%.3f", mort_rate)
    ),
    vjust = 4,       # push slightly above box
    #angle = 45,         # rotate label vertically
    size = 2.5,
    color = "red"
  ) +
  labs(
    title = "Boxplot of Mortality RSI by AHRQ CCS (with Mortality Rate)",
    x = "AHRQ CCS",
    y = "Mortality RSI"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Mort30 and other categorical variables
# Remove outcome variables + ahrq_ccs
plot_cat_vars <- cat_vars[!cat_vars %in% c(vars_outcomes, "ahrq_ccs")]

# Summarize mortality by categorical variables
mort_cat_summary <- data %>%
  # convert binary to factor for plotting
  mutate(across(all_of(vars_comorbidity), ~ as.factor(.))) %>%
  select(all_of(plot_cat_vars), mort30) %>%
  pivot_longer(
    cols      = all_of(plot_cat_vars),
    names_to  = "variable",
    values_to = "level"
  ) %>%
  group_by(variable, level) %>%
  summarise(
    mort_rate = mean(mort30, na.rm = TRUE),
    n         = n(),
    .groups   = "drop"
  )

# Plot
ggplot(mort_cat_summary,
       aes(x = level, y = mort_rate)) +
  geom_col() +
  geom_text(
    aes(label = sprintf("%.3f", mort_rate)),
    vjust = -0.3,
    size  = 3
  ) +
  facet_wrap(~ variable, scales = "free_x") +
  labs(
    title = "30-day mortality rate by categorical variables",
    x     = "Category level",
    y     = "Mortality rate"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

plot_rsi_box_all <- function(data, cat_vars) {
  
  df_long <- data %>%
    mutate(across(all_of(cat_vars), as.factor)) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(cat_vars),
      names_to = "variable",
      values_to = "category"
    )
  
  ggplot(df_long, aes(x = category, y = mortality_rsi)) +
    geom_boxplot() +
    facet_wrap(~ variable, scales = "free_x") +
    labs(
      title = "Mortality RSI by All Categorical Variables",
      x     = "Category",
      y     = "Mortality RSI"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text  = element_text(size = 10)
    )
}

plot_rsi_box_all(data, c("gender", "baseline_cancer", "baseline_dementia","baseline_osteoart","baseline_psych",
                         "baseline_pulmonary"))

############################################################
## relationship between mortality30 and numerical variables
############################################################
# mort30 and mortality_rsi
# higher mortality_rsi should correspond to higher mort30
ggplot(data, aes(x = mortality_rsi, fill = factor(mort30))) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Density plot of Mortality RSI by 30-day mortality",
    x = "Mortality RSI",
    y = "Density",
    fill = "mort30"
  ) +
  theme_minimal()
num_vars <- setdiff(vars_continuous, "mort30")   # remove outcome
df_num_long <- data %>%
  select(all_of(num_vars), mort30) %>%
  pivot_longer(
    cols = all_of(num_vars),
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(mort30 = factor(mort30))
means_df <- df_num_long %>%
  group_by(variable, mort30) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")
ggplot(df_num_long, aes(x = value, fill = mort30)) +
  geom_density(alpha = 0.5) +
  geom_vline(
    data = means_df,
    aes(xintercept = mean_value, color = mort30),
    linetype = "dashed",
    size = 0.7,
    alpha = 0.8
  ) +
  facet_wrap(~ variable, scales = "free") +
  labs(
    title = "Density distributions of numerical variables by 30-day mortality (with mean lines)",
    x = "Value",
    y = "Density",
    fill = "mort30",
    color = "mort30"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10),
    legend.position = "top"
  )


# boxplot and violin plot
ggplot(df_num_long, aes(x = mort30, y = value)) +
  geom_boxplot(outlier.alpha = 0.2) +
  stat_summary(
    fun = mean,
    geom = "point",
    color = "red",
    size = 2,
    position = position_dodge(width = 0.75)
  ) +
  facet_wrap(~ variable, scales = "free") +
  labs(
    title = "Boxplots of Numerical Variables by 30-day Mortality (with Mean Points)",
    x = "Mortality (0 = alive, 1 = dead)",
    y = "Value"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )



# violin plot
ggplot(df_num_long, aes(x = mort30, y = value)) +
  geom_violin(fill = "skyblue", alpha = 0.5) +
  facet_wrap(~ variable, scales = "free") +
  labs(
    title = "Violin Plots of Numerical Variables by 30-day Mortality",
    x = "Mortality (0 = alive, 1 = dead)",
    y = "Value"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# binned
df_num_long2 <- data %>%
  select(all_of(num_vars), mort30) %>%
  pivot_longer(
    cols      = all_of(num_vars),
    names_to  = "variable",
    values_to = "value"
  )
df_num_binned <- df_num_long2 %>%
  group_by(variable) %>%
  mutate(bin = ntile(value, 10)) %>%   # deciles
  group_by(variable, bin) %>%
  summarise(
    mort_rate = mean(mort30, na.rm = TRUE),
    n         = n(),
    value_mid = mean(value, na.rm = TRUE),
    .groups   = "drop"
  )

ggplot(df_num_binned, aes(x = value_mid, y = mort_rate)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ variable, scales = "free_x") +
  labs(
    title = "Binned crude 30-day mortality rate (10 bins) for all numeric variables",
    x = "Binned value (midpoint)",
    y = "Mortality rate"
  ) +
  theme_minimal()

# loess smoothed
ggplot(df_num_long2, aes(x = value, y = mort30)) +
  geom_smooth(method = "loess", se = TRUE) +
  facet_wrap(~ variable, scales = "free") +
  labs(
    title = "Smoothed (LOESS) relationship between numeric variables and 30-day mortality",
    x = "Value",
    y = "Mortality proportion"
  ) +
  theme_minimal()

############################################################
## Distribution of continuous vars by mort30 outcome (facets)
############################################################

# 1. Exclude mortality_rsi (outcome on y-axis)
vars_predictors <- setdiff(vars_continuous, "mortality_rsi")

# 2. Pivot only the predictor variables
df_long <- data %>%
  select(all_of(vars_predictors), mortality_rsi, mort30) %>%
  pivot_longer(
    cols = all_of(vars_predictors),
    names_to = "variable",
    values_to = "value"
  )

# 3. Faceted scatterplot:
ggplot(df_long, aes(x = value, y = mortality_rsi, color = factor(mort30))) +
  geom_point(alpha = 0.45) +
  facet_wrap(~ variable, scales = "free_x") +
  labs(
    title = "Scatterplots of mortality_rsi vs continuous predictors",
    x = "Predictor Value",
    y = "Mortality RSI",
    color = "mort30"
  ) +
  theme_minimal()


############################################################
## 11. Quick summary tables (for paper)
############################################################

# Overall descriptive summary of key variables
data %>%
  summarise(
    n = n(),
    mort30_rate = mean(mort30, na.rm = TRUE),
    comp_rate   = mean(complication, na.rm = TRUE),
    age_mean    = mean(age, na.rm = TRUE),
    age_sd      = sd(age, na.rm = TRUE),
    bmi_mean    = mean(bmi, na.rm = TRUE),
    bmi_sd      = sd(bmi, na.rm = TRUE)
  )

# Mortality by quintiles of hour
data %>%
  mutate(hour_q = ntile(hour, 5)) %>%
  group_by(hour_q) %>%
  summarise(
    hour_min = min(hour, na.rm = TRUE),
    hour_max = max(hour, na.rm = TRUE),
    mort_rate = mean(mort30, na.rm = TRUE),
    n = n()
  )

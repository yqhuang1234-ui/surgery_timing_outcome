# set up the environment
setwd("~/Dropbox/School/CU/fall 2025/BIOS 6618 adv biostatistical method/final project/surgery_timing_outcome")
source("./code/0.0_setup_surgery-timing-outcome.R")

cat("=== Structure ===\n")
str(data)

#############################################################
## association tests between predictors and mort30
#############################################################
#==== Chi-square tests for categorical variables ====#
# exclude outcome variables from categorical variable list
vars_categorical_binary <- vars_categorical_binary[!vars_categorical_binary %in% vars_outcomes]
# run chi-square tests for categorical variables
chi_square_results <- data.frame(
  variable   = character(),
  category   = character(),
  mort_rate  = numeric(),
  p_value    = numeric(),
  sig_level  = character(),
  stringsAsFactors = FALSE
)
for (v in vars_categorical_binary) {
  tab <- table(data[[v]], data$mort30)
  chi_test <- chisq.test(tab)
  p_val <- chi_test$p.value
  mort_rates <- prop.table(tab, 1)[, "1"]
  
  for (cat in rownames(tab)) {
    chi_square_results <- rbind(
      chi_square_results,
      data.frame(
        variable  = v,
        category  = cat,
        mort_rate = mort_rates[cat],
        p_value   = p_val,
        sig_level = ifelse(p_val < 0.05, "Yes", "No"),
        stringsAsFactors = FALSE
      )
    )
  }
}
# pivot the table showing category value as columns for better readability
chi_square_results <- chi_square_results %>%
  tidyr::pivot_wider(
    names_from  = category,
    values_from = mort_rate,
    names_prefix = "mort_rate_"
  )
print(chi_square_results)
# return variables with significant association
sig_cat_vars <- unique(chi_square_results$variable[chi_square_results$sig_level == "Yes"])
cat("Variables with significant association with mort30 (p < 0.05):\n")
print(sig_cat_vars)

#====t test for continuous variables ====#
t_test_results <- data.frame(
  variable   = character(),
  mean_mort0 = numeric(),
  mean_mort1 = numeric(),
  p_value    = numeric(),
  sig_level  = character(),
  stringsAsFactors = FALSE
)
for (v in vars_continuous) {
  group0 <- data[[v]][data$mort30 == 0]
  group1 <- data[[v]][data$mort30 == 1]
  t_test <- t.test(group0, group1)
  p_val <- t_test$p.value
  mean0 <- mean(group0, na.rm = TRUE)
  mean1 <- mean(group1, na.rm = TRUE)
  t_test_results <- rbind(
    t_test_results,
    data.frame(
      variable   = v,
      mean_mort0 = mean0,
      mean_mort1 = mean1,
      p_value    = p_val,
      sig_level  = ifelse(p_val < 0.05, "Yes", "No"),
      stringsAsFactors = FALSE
    )
  )
}
print(t_test_results)
# return variables with significant association
sig_cont_vars <- t_test_results$variable[t_test_results$sig_level == "Yes"]
cat("Continuous variables with significant difference by mort30 (p < 0.05):\n")
print(sig_cont_vars)

# conduct wilcoxon test for continuous variables that are not normally distributed
# check normality using shapiro-wilk test skip variables with error and include in the list of non-normal variables
non_normal_vars <- c()
for (v in vars_continuous) {
  group0 <- data[[v]][data$mort30 == 0]
  group1 <- data[[v]][data$mort30 == 1]
  shapiro0 <- tryCatch(shapiro.test(group0), error = function(e) NULL)
  shapiro1 <- tryCatch(shapiro.test(group1), error = function(e) NULL)
  if (!is.null(shapiro0) && shapiro0$p.value < 0.05) {
    non_normal_vars <- c(non_normal_vars, v)
    next
  }
  if (!is.null(shapiro1) && shapiro1$p.value < 0.05) {
    non_normal_vars <- c(non_normal_vars, v)
    next
  }
}
cat("Continuous variables that are not normally distributed (p < 0.05):\n")
print(non_normal_vars)
# run wilcoxon test for non-normal continuous variables
wilcoxon_test_results <- data.frame(
  variable   = character(),
  median_mort0 = numeric(),
  median_mort1 = numeric(),
  p_value    = numeric(),
  sig_level  = character(),
  stringsAsFactors = FALSE
)
for (v in non_normal_vars) {
  group0 <- data[[v]][data$mort30 == 0]
  group1 <- data[[v]][data$mort30 == 1]
  wilcoxon_test <- wilcox.test(group0, group1)
  p_val <- wilcoxon_test$p.value
  median0 <- median(group0, na.rm = TRUE)
  median1 <- median(group1, na.rm = TRUE)
  wilcoxon_test_results <- rbind(
    wilcoxon_test_results,
    data.frame(
      variable     = v,
      median_mort0 = median0,
      median_mort1 = median1,
      p_value      = p_val,
      sig_level    = ifelse(p_val < 0.05, "Yes", "No"),
      stringsAsFactors = FALSE
    )
  )
}
print(wilcoxon_test_results)
# return variables with significant association
sig_wilcoxon_vars <- wilcoxon_test_results$variable[wilcoxon_test_results$sig_level == "Yes"]
cat("Continuous variables with significant difference by mort30 (p < 0.05) from Wilcoxon test:\n")
print(sig_wilcoxon_vars)
print(t_test_results)


# auto-select test

auto_test_continuous <- function(data, outcome, vars) {
  out <- list()
  
  # check outcome exists
  if (!outcome %in% names(data)) {
    stop("Outcome variable not found in data.")
  }
  
  for (v in vars) {
    if (!v %in% names(data)) {
      warning(paste("Variable", v, "not found in data. Skipping."))
      next
    }
    
    x <- data[[v]]
    y <- data[[outcome]]
    
    # remove rows with NA in x or y
    complete_idx <- complete.cases(x, y)
    x <- x[complete_idx]
    y <- y[complete_idx]
    
    # split into groups: assume outcome is 0/1
    g1 <- x[y == 0]
    g2 <- x[y == 1]
    
    n1 <- length(g1)
    n2 <- length(g2)
    
    # default
    method <- NA
    pval   <- NA
    
    if (n1 < 2 || n2 < 2) {
      # not enough data to compare
      method <- "insufficient data"
      pval   <- NA
    } else if (n1 >= 30 && n2 >= 30) {
      # large samples: t-test OK by CLT
      method <- "t-test (large N)"
      test_res <- t.test(g1, g2)  # Welch t-test by default
      pval <- test_res$p.value
    } else {
      # small samples: try Shapiro when size allows
      p_norm1 <- if (n1 >= 3 && n1 <= 5000) shapiro.test(g1)$p.value else NA
      p_norm2 <- if (n2 >= 3 && n2 <= 5000) shapiro.test(g2)$p.value else NA
      
      if (!is.na(p_norm1) && !is.na(p_norm2) &&
          p_norm1 > 0.05 && p_norm2 > 0.05) {
        method <- "t-test (normality OK)"
        test_res <- t.test(g1, g2)
        pval <- test_res$p.value
      } else {
        method <- "Wilcoxon (Mann-Whitney U)"
        test_res <- wilcox.test(g1, g2, exact = FALSE)
        pval <- test_res$p.value
      }
    }
    
    out[[v]] <- data.frame(
      variable  = v,
      method    = method,
      p_value   = pval,
      mean_0    = mean(g1, na.rm = TRUE),
      mean_1    = mean(g2, na.rm = TRUE),
      median_0  = median(g1, na.rm = TRUE),
      median_1  = median(g2, na.rm = TRUE),
      n0        = n1,
      n1        = n2,
      stringsAsFactors = FALSE
    )
  }
  
  bind_rows(out)
}

results <- auto_test_continuous(data, outcome = "mort30", vars = vars_continuous)
print(results)


#############################################################
## multicollinearity check for continuous variables
#############################################################
# calculate vif for continuous variables
vif_data <- data %>%
  select(all_of(c(vars_continuous, "mort30"))) %>%
  na.omit()
vif_model <- lm(mort30 ~ ., data = vif_data)
vif_values <- vif(vif_model)
vif_results <- data.frame(
  variable = names(vif_values),
  vif      = as.numeric(vif_values),
  stringsAsFactors = FALSE
)
print(vif_results)
# return variables with vif > 5
high_vif_vars <- vif_results$variable[vif_results$vif > 10]
cat("Continuous variables with high VIF (> 5):\n")
print(high_vif_vars)

# calculate correlation matrix for continuous variables
cor_data <- data %>%
  select(all_of(vars_continuous)) %>%
  na.omit()
cor_matrix <- cor(cor_data)
print(cor_matrix)
# return variable pairs with high correlation (> 0.7)
# remove duplicate pairs and self-correlation
high_cor_results <- data.frame(
  variable1 = character(),
  variable2 = character(),
  correlation = numeric(),
  stringsAsFactors = FALSE
)
for (i in 1:(ncol(cor_matrix) - 1)) {
  for (j in (i + 1):ncol(cor_matrix)) {
    cor_value <- cor_matrix[i, j]
    if (abs(cor_value) > 0.7) {
      high_cor_results <- rbind(
        high_cor_results,
        data.frame(
          variable1   = colnames(cor_matrix)[i],
          variable2   = colnames(cor_matrix)[j],
          correlation = cor_value,
          stringsAsFactors = FALSE
        )
      )
    }
  }
}

# Based on the results, consider removing one variable from each highly correlated pair
numeric_vars_remove <- c("bmi","complication_rsi","ccsComplicationRate")
# final continuous variables to consider
numeric_vars_consider <- setdiff(vars_continuous, numeric_vars_remove)
cat("Variable pairs with high correlation (> 0.7):\n")
print(high_cor_results)
cat("Removing the following continuous variables due to high multicollinearity:\n")
print(numeric_vars_remove)
cat("Continuous variables to consider after removing multicollinear variables:\n")
print(numeric_vars_consider)

#############################################################
#missing value imputation
#############################################################
# check missingness
missing_summary <- data.frame(
  variable       = character(),
  missing_count  = numeric(),
  missing_percent = numeric(),
  stringsAsFactors = FALSE
)
for (v in colnames(data)) {
  missing_count <- sum(is.na(data[[v]]))
  missing_percent <- missing_count / nrow(data) * 100
  missing_summary <- rbind(
    missing_summary,
    data.frame(
      variable        = v,
      missing_count   = missing_count,
      missing_percent = missing_percent,
      stringsAsFactors = FALSE
    )
  )
}
# check missingness by mort30 and pivot by mort30
missing_by_mort30 <- data %>%
  group_by(mort30) %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  tidyr::pivot_longer(
    cols = -mort30,
    names_to = "variable",
    values_to = "missing_count"
  ) %>%
  tidyr::pivot_wider(
    names_from = mort30,
    values_from = missing_count,
    names_prefix = "mort30_"
  )
missing_summary <- missing_summary %>%
  left_join(missing_by_mort30, by = "variable")
print(missing_summary)

#############################################
# finalize variable lists for modeling
#############################################
vars_to_check <- c("age", "gender", "asa_status")
vars_to_drop  <- c("bmi", 
                   "ccsComplicationRate","complication_rsi",
                   vars_categorical_binary[!vars_categorical_binary %in% sig_cat_vars]) 
cat("Variables to drop due to non-significant association with mort30 or high correlation or high missingness:\n")
print(vars_to_drop)
data_clean <- data %>%
  drop_na(all_of(vars_to_check)) %>%
  select(-all_of(vars_to_drop))
cat("Data dimensions after dropping missing values in key variables and unnecessary variables:\n")
print(dim(data_clean))
print(length(vars_to_drop))
length(vars_to_drop) + ncol(data_clean) == ncol(data)  # should equal ncol(data)
# final variable lists
final_vars_categorical <- setdiff(vars_categorical, vars_to_drop)
final_vars_binary      <- setdiff(vars_binary, vars_to_drop)
final_vars_continuous  <- numeric_vars_consider
cat("Final categorical variables:\n")
print(final_vars_categorical)
cat("Final binary variables:\n")
print(final_vars_binary)
cat("Final continuous variables:\n")
print(final_vars_continuous)
# check missingness again
final_missing_summary <- data.frame(
  variable       = character(),
  missing_count  = numeric(),
  missing_percent = numeric(),
  stringsAsFactors = FALSE
)
for (v in colnames(data_clean)) {
  missing_count <- sum(is.na(data_clean[[v]]))
  missing_percent <- missing_count / nrow(data_clean) * 100
  final_missing_summary <- rbind(
    final_missing_summary,
    data.frame(
      variable        = v,
      missing_count   = missing_count,
      missing_percent = missing_percent,
      stringsAsFactors = FALSE
    )
  )
}
print(final_missing_summary)

#############################################################
# combine ahrq_ccs levels for modeling
#############################################################
data_clean <- data_clean %>%
  mutate(
    procedure_group = case_when(
      ahrq_ccs %in% c("Colorectal resection",
                      "Gastrectomy; partial and total",
                      "Small bowel resection",
                      "Inguinal and femoral hernia repair",
                      "Other hernia repair") ~ "Gastrointestinal & Hernia",
      
      ahrq_ccs %in% c("Lumpectomy; quadrantectomy of breast",
                      "Mastectomy") ~ "Breast",
      
      ahrq_ccs %in% c("Arthroplasty knee",
                      "Hip replacement; total and partial",
                      "Laminectomy; excision intervertebral disc",
                      "Spinal fusion") ~ "Orthopedic & Spine",
      
      ahrq_ccs %in% c("Hysterectomy; abdominal and vaginal",
                      "Oophorectomy; unilateral and bilateral",
                      "Other excision of cervix and uterus",
                      "Repair of cystocele and rectocele; obliteration of vaginal vault") ~ "Gynecologic",
      
      ahrq_ccs %in% c("Endoscopy and endoscopic biopsy of the urinary tract",
                      "Nephrectomy; partial or complete",
                      "Open prostatectomy",
                      "Transurethral resection of prostate (TURP)") ~ "Urologic",
      
      ahrq_ccs %in% c("Thyroidectomy; partial or complete") ~ "Endocrine",
      
      ahrq_ccs %in% c("Plastic procedures on nose") ~ "Plastic",
      
      TRUE ~ "Other non-major procedures"
    )
  )
# check procedure group by original ahrq_ccs
print(unique(data_clean$procedure_group))
procedure_check <- data_clean %>%
  group_by(procedure_group,ahrq_ccs) %>%
  summarise(count = n()) %>%
  arrange(procedure_group)
print(procedure_check, n=length(unique(data_clean$ahrq_ccs)))

#plot histogram of mort30 by procedure_group
ggplot(data_clean, aes(x=procedure_group, fill=factor(mort30))) +
  geom_bar(position="dodge") +
  labs(title="30-day Mortality by Procedure Group", x="Procedure Group", y="Count", fill="30-day Mortality") +
  theme_minimal()

# summarize mort30 by procedure_group
# check sample size (n<100?) and death(n<10?) in each procedure_group
group_summary <- data_clean %>%
  group_by(procedure_group) %>%
  summarise(
    n          = n(),
    deaths     = sum(mort30, na.rm = TRUE),
    mort_rate  = mean(mort30, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  arrange(desc(n))

group_summary

# check heterogeneity within each group
# within the group, are mortality rates similar across different ahrq_ccs?
within_group <- data_clean %>%
  group_by(procedure_group, ahrq_ccs) %>%
  summarise(
    n          = n(),
    deaths     = sum(mort30, na.rm = TRUE),
    mort_rate  = mean(mort30, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  arrange(procedure_group, desc(n))

print(within_group, n=27)

ggplot(group_summary, aes(x = procedure_group, y = mort_rate)) +
  geom_col() +
  geom_text(
    aes(label = sprintf("%.2f%%", mort_rate * 100)),
    vjust = -0.3,
    size  = 3
  ) +
  labs(
    title = "30-day mortality rate by procedure group",
    x     = "Procedure group",
    y     = "Mortality rate"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# combine group with < 6 deaths
data_clean <- data_clean %>%
  mutate(
    procedure = ifelse(
      procedure_group %in% group_summary$procedure_group[group_summary$deaths < 6],
      "Other non-major procedures",
      procedure_group
    )
  )
# check new procedure group summary
new_group_summary <- data_clean %>%
  group_by(procedure) %>%
  summarise(
    n          = n(),
    deaths     = sum(mort30, na.rm = TRUE),
    mort_rate  = mean(mort30, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  arrange(desc(n))
print(new_group_summary)

#check finalized group, procedure group and ahrq_ccs
final_procedure_check <- data_clean %>%
  group_by(procedure, procedure_group,ahrq_ccs) %>%
  summarise(count = n()) %>%
  arrange(procedure)
print(final_procedure_check, n=length(unique(data_clean$ahrq_ccs)))

####################################################
#encode final categorical variables
####################################################
# replace ahrq_ccs with procedure_group2
final_vars_categorical <- setdiff(final_vars_categorical, "ahrq_ccs")
final_vars_categorical <- c(final_vars_categorical, "procedure")
print(final_vars_categorical)
# check unique levels in each categorical variable
for (v in final_vars_categorical) {
  cat("Variable:", v, "\n")
  print(unique(data_clean[[v]]))
}
data_clean <- data_clean %>%
  mutate(
    procedure = factor(
      procedure,
      levels = c(
        "Other non-major procedures",  # reference
        "Gynecologic",
        "Urologic",
        "Orthopedic & Spine",
        "Gastrointestinal & Hernia"
      )
    )
  )
# drop original procedure_group and ahrq_ccs
data_to_model <- data_clean %>%
  select(-c(procedure_group, ahrq_ccs))
# encode categorical variables using one-hot encoding
data_model <- data_to_model%>%
  mutate(across(all_of(final_vars_categorical), as.factor)) %>%
  fastDummies::dummy_cols(
    select_columns = final_vars_categorical,
    remove_first_dummy = TRUE,
    remove_selected_columns = TRUE
  )
cat("Data dimensions after one-hot encoding categorical variables:\n")
print(dim(data_model))
# check new variable names
new_categorical_vars <- colnames(data_model)[grepl(paste0(final_vars_categorical, collapse = "|"), colnames(data_model))]
cat("New categorical variable names after one-hot encoding:\n")
print(new_categorical_vars)
names(data_model)

# not encode for asa status since it is ordinal
# encode categorical variables using one-hot encoding
data_model_no_asa <- data_to_model %>%
  mutate(across(all_of(c("gender",
                         "procedure")), as.factor)) %>%
  fastDummies::dummy_cols(
    select_columns = c("gender",
                       "procedure"),
    remove_first_dummy = TRUE,
    remove_selected_columns = TRUE
  )
cat("Data dimensions after one-hot encoding categorical variables:\n")
print(dim(data_model))
# check new variable names
new_categorical_vars_no_asa <- colnames(data_model_no_asa)[grepl(paste0(c("gender","procedure"), collapse = "|"), colnames(data_model_no_asa))]
cat("New categorical variable names after one-hot encoding:\n")
print(new_categorical_vars_no_asa)
names(data_model_no_asa)
str(data_model_no_asa)

str(data_model)
###################################################
# save cleaned data for modeling
###################################################
#add today's date to file name
today_date <- format(Sys.Date(), "%Y-%m-%d")
print(today_date)
#create directory if not exist
if (!dir.exists("./data/processed/")) {
  dir.create("./data/processed/")
}
#save data
file_name_asa_factor <- paste0(today_date, "_dropna_regroup-procedure_encode-factors_asa_factor.rds")
saveRDS(data_model, file = paste0("./data/processed/", file_name_asa_factor))
file_name_asa_ordinal <- paste0(today_date, "_dropna_regroup-procedure_encode-factors_asa_ordinal.rds")
saveRDS(data_model_no_asa, file = paste0("./data/processed/", file_name_asa_ordinal))
filename_without_encode <- paste0(today_date, "_dropna_regroup-procedure_no-encode.rds")
saveRDS(data_to_model, file = paste0("./data/processed/", filename_without_encode))
cat("Cleaned data saved for modeling:\n")
print(paste0("./data/processed/", file_name_asa_factor))
print(paste0("./data/processed/", file_name_asa_ordinal))
print(paste0("./data/processed/", filename_without_encode))

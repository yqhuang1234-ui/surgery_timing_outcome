# set up the environment
setwd("~/Dropbox/School/CU/fall 2025/BIOS 6618 adv biostatistical method/final project/surgery_timing_outcome")
source("./code/0.0_setup_surgery-timing-outcome.R")

cat("=== Structure ===\n")
str(data)


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
print(high_cor_results)
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
vars_to_drop  <- c("bmi", "race", 
                   "ccsComplicationRate","complication_rsi"
                   ) 
cat("Variables to drop due to  high correlation or high missingness:\n")
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
    ),
    dow = factor(
      dow,
      levels=c("1","2","3","4","5")
    )
  )
# drop original procedure_group and ahrq_ccs
data_to_model <- data_clean %>%
  select(-c(procedure_group, ahrq_ccs, "complication","moonphase"))
str(data_to_model)

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
filename_without_encode <- paste0(today_date, "_dropna_regroup-procedure_no-encode.rds")
saveRDS(data_to_model, file = paste0("./data/processed/", filename_without_encode))
cat("Cleaned data saved for modeling:\n")
print(paste0("./data/processed/", filename_without_encode))

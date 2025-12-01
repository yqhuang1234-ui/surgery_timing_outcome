# set up the environment
setwd("~/Dropbox/School/CU/fall 2025/BIOS 6618 adv biostatistical method/final project/surgery_timing_outcome")
source("./code/0.0_setup_surgery-timing-outcome.R")

cat("=== Structure ===\n")
str(data)
summary(data)

#############################################################
# logic cleaning
#############################################################
# keep age >= 18
orin_dim <- dim(data)
# find mort30 for age < 18
mort_young <- data %>%
  filter(age < 18) %>%
  select(mort30)
cat("30-day mortality for patients age < 18:\n")
print(table(mort_young$mort30))
print(paste0("Original data dimensions: ", orin_dim[1], " rows, ", orin_dim[2], " columns"))
data <- data %>%
  filter(age >= 18)
cat("Data dimensions after filtering age >= 18:\n")
print(dim(data))
removed_rows <- orin_dim[1] - nrow(data)
cat(paste0("Number of rows removed due to age < 18: ", removed_rows, "\n"))
# check age summary
summary(data$age)


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
    ),
    month = factor(
      month,
      levels=c("1","2","3","4","5","6","7","8","9","10","11","12")
    )
  )
# drop original procedure_group and ahrq_ccs
data_to_model <- data_clean %>%
  select(-c(procedure_group, ahrq_ccs, "complication","moonphase","dow","month","ccsMort30Rate"))
str(data_to_model)

# -------------------------
# 1. CREATE HOUR CATEGORIES
# -------------------------
data_to_model$hour_cat <- cut(
  data_to_model$hour,
  breaks = c(6, 8, 11, 13, 15, 17, 19),
  right = FALSE,  # intervals left-inclusive
  include.lowest = FALSE,
  labels = c("6-8am", "8-11am", "11-1pm", "1-3pm", "3-5pm", "5-7pm")
)

# check distribution
table(data_to_model$hour_cat, useNA = "ifany")
#return row for na hour_cat
data_to_model[is.na(data_to_model$hour_cat), ]
#if hour=19 set it to "5-7"
data_to_model$hour_cat[is.na(data_to_model$hour_cat) & data_to_model$hour == 19] <- "5-7pm"
table(data_to_model$hour_cat, useNA = "ifany")



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


# baseline characteristics summaries between hour group
# create figure directory if not exist
figure_dir <- file.path("results", today_date,'/')
print(figure_dir)

if (!dir.exists(figure_dir)) {
  dir.create(figure_dir, recursive = TRUE)
}
## --------------------------------------------------
## 0. Pre-processing: factors, labels, 0/1 → No/Yes
## --------------------------------------------------

# set your hour_cat levels if you haven't already
data_to_plot <- data_to_model %>%
  mutate(
    # gender: 0 = Male, 1 = Female
    gender = factor(ifelse(gender == 1, "Female", "Male"),
                    levels = c("Male", "Female")),
    mort30 = factor(ifelse(mort30 == 1, "Yes", "No"),
                    levels = c("No", "Yes")),
    asa_status = factor(asa_status,
                        levels = c("1", "2", "3"),
                        labels = c(
                          "I-II",
                          "III",
                          "IV-VI"
                        )
    )
  )

# list of binary baseline comorbidity variables that are 0/1
comorb_vars <- c(
  "baseline_cancer", "baseline_cvd", "baseline_dementia",
  "baseline_diabetes", "baseline_digestive",
  "baseline_osteoart", "baseline_psych", "baseline_pulmonary"
)

# recode all of them 0/1 -> No/Yes
data_to_plot <- data_to_plot %>%
  mutate(
    across(all_of(comorb_vars),
           ~ factor(ifelse(. == 1, "Yes", "No"),
                    levels = c("No", "Yes")))
  )

## --------------------------------------------------
## 1. Variable sets & ordering
## --------------------------------------------------

cont_vars <- c("age", "mortality_rsi", "baseline_charlson")

cat_vars <- names(data_to_plot)[
  !(names(data_to_plot) %in% c(cont_vars,  "hour", "hour_cat"))
]

# put main variables first, then comorbidities, then any other cats
main_cat <- c("gender", "asa_status","procedure")
other_cat <- setdiff(cat_vars, c(main_cat, comorb_vars))

var_order <- c(
  "age",
  "mortality_rsi",
  "baseline_charlson",
  main_cat,
  comorb_vars,
  other_cat
)
var_order
## --------------------------------------------------
## 2. gtsummary baseline table
## --------------------------------------------------

tbl_baseline <-
  data_to_plot %>%
  select(hour_cat, all_of(var_order)) %>%
  tbl_summary(
    by = hour_cat,
    statistic = list(
      all_continuous()  ~ "{mean} ± {sd}",
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "no",
    label = list(
      age              ~ "Age",
      mortality_rsi    ~ "Risk Stratification Index (30-day mortality)",
      baseline_charlson~ "Charlson Comorbidity Index",
      gender           ~ "Gender",
      procedure        ~ "Procedure",
      asa_status       ~ "ASA Physical Status",
      # comorbidities
      baseline_cancer      ~ "Cancer",
      baseline_cvd         ~ "Cardiovascular/Cerebrovascular disease",
      baseline_dementia    ~ "Dementia",
      baseline_diabetes    ~ "Diabetes",
      baseline_digestive   ~ "Digestive disease",
      baseline_osteoart    ~ "Osteoarthritis",
      baseline_psych       ~ "Psychiatric disorder",
      baseline_pulmonary   ~ "Pulmonary disease",
      mort30               ~ "30-day Mortality"
    )
  ) %>%
  add_overall(last = FALSE) %>%   
  add_p(
    test = list(
      all_continuous()  ~ "oneway.test",
      all_categorical() ~ "chisq.test",
      baseline_dementia ~ "fisher.test",
      mort30            ~ "fisher.test"
    ),
    test.args = list(
      all_continuous()  ~ list(var.equal = TRUE),
      baseline_dementia ~ list(simulate.p.value = TRUE, B = 10000),
      mort30            ~ list(simulate.p.value = TRUE, B = 10000)
    )
  ) %>%
  modify_header(label ~ "**Variable**") %>%
  bold_labels()

## --------------------------------------------------
## 3. Convert to gt and save PNG
##    (wide but not too tall)
## --------------------------------------------------

gt_tbl <-
  as_gt(tbl_baseline) %>%
  gt::tab_options(
    table.font.size   = 9,           # smaller font
    data_row.padding  = px(1),       # tighten row spacing
    heading.padding   = px(2)
  )
gt_tbl

gtsave(
  gt_tbl,
  filename = paste0(figure_dir, "1_baseline-hour-group-summary_table.png"),
  vwidth  = 2400,  # pixels: wide
  vheight = 700,   # pixels: relatively short
  zoom=6
)

# extract p-values table
sig_vars <- tbl_baseline$table_body %>%
  select(variable, p.value) %>%        # these columns DO exist in table_body
  filter(!is.na(p.value), p.value < 0.05)

sig_vars

# Convert continuous to 4 quantile groups for meaningful comparison
#rows
(32001-212-2-8) == nrow(data_to_model)
names(data_to_model)
get_uni_prop <- function(df, var) {
  var_sym <- rlang::ensym(var)
  
  df %>%
    group_by(!!var_sym) %>%
    summarise(
      mort_rate = mean(mort30),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(
      variable = rlang::as_string(var_sym),
      category = as.character(!!var_sym)
    )
}

data_uni <- data_to_model %>%
  mutate(
    age = cut2(age, g = 4),
    mortality_rsi = cut2(mortality_rsi, g = 4),
    hour = cut2(hour, g = 4),
    baseline_charlson = cut2(baseline_charlson, g = 4)
  )

# filter vars
remove_vec <- c("mort30","hour")
vars <- setdiff(names(data_to_model), remove_vec)

uni_list <- lapply(vars, function(v) get_uni_prop(data_uni, !!sym(v)))
uni <- bind_rows(uni_list)

main_vars <- c("age", "hour_cat", "asa_status", "mortality_rsi",
               "baseline_charlson", "procedure", "gender")

uni_main <- uni %>% 
  filter(variable %in% main_vars) %>%
  mutate(variable = factor(variable, levels = main_vars))

uni_base <- uni %>% 
  filter(!variable %in% main_vars) %>%
  mutate(variable = factor(variable, levels = sort(unique(variable))))
xmax <- max(uni$mort_rate, na.rm = TRUE)
p_main <- ggplot(uni_main, aes(x = mort_rate, y = fct_rev(category))) +
  geom_point(size = 2, alpha = 0.9) +
  facet_wrap(~ variable, scales = "free_y", ncol = 1, strip.position = "top") +
  scale_x_continuous(limits = c(0, xmax)) + 
  labs(
    x = "30-day Mortality Proportion",
    y = NULL,
    title = "Univariable Summaries of 30-day Mortality"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.5, "lines"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
  )
p_base <- ggplot(uni_base, aes(x = mort_rate, y = fct_rev(category))) +
  geom_point(size = 2, alpha = 0.9) +
  facet_wrap(~ variable, scales = "free_y", ncol = 1, strip.position = "top") +
  scale_x_continuous(limits = c(0, xmax)) + 
  labs(
    x = "30-day Mortality Proportion",
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.5, "lines"),
    plot.title = element_blank()
  )
p_main + p_base + plot_layout(ncol = 2, widths = c(1, 1))


# save figure
ggsave(paste0(figure_dir, "1_univariable-summary_figure.png"), width=6, height=8, dpi = 300)
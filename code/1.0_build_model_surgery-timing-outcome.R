#################################
#load data
#################################
# set up the environment
setwd("~/Dropbox/School/CU/fall 2025/BIOS 6618 adv biostatistical method/final project/surgery_timing_outcome")
source("./code/0.0_setup_surgery-timing-outcome.R")
# load processed data
file_name_asa_factor <- "./data/processed/2025-11-24_dropna_regroup-procedure_encode-factors_asa_factor.rds"
data <- readRDS(file_name_asa_factor)
str(data)
# check missingness for data
miss_summary <- naniar::miss_var_summary(data)
print(miss_summary)
names(data)

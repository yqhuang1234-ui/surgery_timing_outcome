## ================================
## 0. Load needed packages and data
## ================================
# Helper function: install if not already available
load_or_install <- function(pkg){
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# List of needed packages
packages <- c("dplyr", "ggplot2",  "kableExtra", "webshot2","magick","knitr","tidyverse","tidyr","ggh4x","GGally")
# Load or install each
lapply(packages, load_or_install)
setwd("~/Library/CloudStorage/Dropbox/School/CU/fall 2025/BIOS 6618 adv biostatistical method/final project/surgery_timing_outcome/")
data <- read.csv("./data/Surgery_Timing.csv")

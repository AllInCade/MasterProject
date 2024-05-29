# List of packages required for the project (may not be 100% complete and/or contain redundancies)
packages <- c("ggplot2", "dplyr", "tidyr", "glmmTMB", "mgcv", "glmnet", "DHARMa", 
              "future","furrr","parallel","zoo","forecast","purrr","gamm4","progressr",
              "TTR","MASS","caret","spdep","sp","gstat","ROSE","pROC","pheatmap" ,"quantmod",
              "future.apply", "haven", "lubridate", "reshape2", "readr", "purrr", "VGAM", "optimx")

# Function to check and install missing packages
install.packages.if.missing <- function(packages) {
  missing_packages <- packages[!packages %in% installed.packages()[,"Package"]]
  if(length(missing_packages)) {
    install.packages(missing_packages, dependencies = TRUE)
  }
}

# Install missing packages
install.packages.if.missing(packages)

# Load packages
lapply(packages, library, character.only = TRUE)

# List of packages required for the project
packages <- c("ggplot2", "dplyr", "tidyr", "glmmTMB", "mgcv")

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

## Installation of Required Packages

source("requirements.R")

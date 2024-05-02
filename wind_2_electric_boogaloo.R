#remotes::install_github("https://github.com/glmmTMB/glmmTMB/tree/smooth_fixes")
library(parallel)
library(doParallel)
library(progress)
library(dplyr)
library(zoo)
library(lubridate)
library(readr)
library(mgcv)
library(glmmTMB)
library(randomForest)
library(forecast)
library(ggplot2)
library(DHARMa)
library(reshape2)
library(purrr)
library(furrr)
library(future)
library(future.apply)
library(gamm4)
library(caret)
library(pROC)
library(precrec)
library(pheatmap)
library(CalibratR)
library(gains)
library(spdep)
library(sp)
library(gstat)
library(ROSE)
library(tidyverse)
library(lubridate)

plan(sequential)
gc()
options(future.seed = TRUE, future.rng.onMisuse="ignore")

# Read the data
wind_data_full <- read_csv("SB_Aeolus_20230519.csv") %>%
  na.omit() %>%
  mutate(datetime_UTC = ymd_hms(datetime_UTC))

# Filter data for January 2021 and from all stations without NAs
jan_data <- wind_data_full %>%
  filter(year(datetime_UTC) == 2020, month(datetime_UTC) == 1) %>%
  na.omit()

# Prepare lagged data for station 46078
lagged_data <- jan_data %>%
  filter(abb == 46078) %>%
  arrange(datetime_UTC) %>%
  mutate(
    time = row_number(),
    spd_lag = lag(spd),
    U_lag = lag(U),
    V_lag = lag(V),
    U10_lag = lag(U10),
    V10_lag = lag(V10),
    Utau_lag = lag(Utau),
    Vtau_lag = lag(Vtau),
    Atmp_lag = lag(Atmp),
    Wtmp_lag = lag(Wtmp),
    abar_lag = lag(abar),
    dir_lag = lag(dir)
  ) %>%
  filter(!is.na(time))

# Real-time data from other stations
real_time_data <- jan_data %>%
  filter(abb != 46078)

# Combining lagged data from station 46078 with real-time data from other stations
combined_data <- lagged_data %>%
  left_join(real_time_data, by = "datetime_UTC", suffix = c("_lag", "_rt"))


rm(wind_data_full)
gc()
wind_data <- wind_data[order(wind_data$time), ]
#n <- nrow(wind_data)

#wind_ts <- ts(wind_data$spd, frequency = 24) # capture diurnal pattern
#sarima_model <- auto.arima(wind_ts)
#residuals_data <- residuals(sarima_model)
#wind_data$sarima_res <- residuals_data

#wind_data <- wind_data %>%
#  arrange(datetime_UTC) %>% # Make sure the data is sorted by datetime_UTC
#  mutate(time = as.numeric(difftime(datetime_UTC, min(datetime_UTC), units = "hours")) + 1)

wind_data$datetime_UTC <- NULL # date type variables sometimes cause problems with predict()

wind_data <- wind_data %>%
  mutate(
    time = lag(time, 1),
    spd_lag = lag(spd, 1),
    U = lag(U, 1),
    V = lag(V, 1),
    U10 = lag(U10, 1),
    V10 = lag(V10, 1),
    Utau = lag(Utau, 1),
    Vtau = lag(Vtau, 1),
    Atmp = lag(Atmp, 1),
    Wtmp = lag(Wtmp, 1),
    abar = lag(abar, 1),
    dir = lag(dir, 1)
  ) %>%
  filter(!is.na(time)) 

rf_model <- randomForest(spd ~ ., data = wind_data, ntree = 1000, importance = TRUE)
importance_matrix <- importance(rf_model)
sorted_importance_mse <- sort(importance_matrix[, "%IncMSE"], decreasing = TRUE)
sorted_importance_purity <- sort(importance_matrix[, "IncNodePurity"], decreasing = TRUE)
print(sorted_importance_mse)
print(sorted_importance_purity)

numeric_df <- wind_data[sapply(wind_data, is.numeric)]
x <- as.matrix(numeric_df[, setdiff(names(numeric_df), "spd")])
y <- numeric_df$spd
set.seed(123)  # for reproducibility
lasso_model <- glmnet(x, y, alpha = 1)  # alpha = 1 for lasso
plot(lasso_model, xvar = "lambda", label = TRUE)
cv_model <- cv.glmnet(x, y, alpha = 1)
plot(cv_model)
best_lambda <- cv_model$lambda.min
best_model <- glmnet(x, y, alpha = 1, lambda = best_lambda)
coef(best_model)

calculate_rmse <- function(actual, predicted) {
  sqrt(mean((actual - predicted)^2))
}

plan(sequential)
gc()
options(future.seed = TRUE, future.rng.onMisuse="ignore")

total_obs <- nrow(wind_data)
training_window_size <- floor(0.9*total_obs)
forecast_horizon <- 24 
plan(multisession, workers = 6) # potential RAM bottleneck

process_fold <- function(start_index, data, training_window_size, forecast_horizon) {
  train_end <- start_index + training_window_size - 1
  test_end <- train_end + forecast_horizon
  
  train_data <- data[start_index:train_end, ]
  test_data <- data[(train_end + 1):test_end, ]
  
  lin_model <- lm(spd ~ time + spd_lag + abar + dir + Atmp + Wtmp,
                        data = train_data)
  
  gam_model <- glmmTMB(spd ~ s(time) + s(spd_lag) + 
                       abar + Atmp + Wtmp + dir,
                       disp =~ Utau,
                       data = train_data,
                       family = gaussian(), 
                       REML = TRUE)
  
  lin_predictions <- predict(lin_model, newdata = test_data)
  gam_predictions <- predict(gam_model, newdata = test_data, type = "response")
  
  lin_rmse <- calculate_rmse(test_data$spd, lin_predictions)
  gam_rmse <- calculate_rmse(test_data$spd, gam_predictions)
  
  return(list(
    lin_rmse = lin_rmse,
    gam_rmse = gam_rmse,
    actual_spd = test_data$spd,
    lin_predictions = lin_predictions,
    gam_predictions = gam_predictions,
    fold = rep(start_index, nrow(test_data))  # Use 'start_index' instead of 'fold_index'
  ))
}

start_indices <- 1:(nrow(wind_data) - training_window_size - forecast_horizon + 1)
results <- future_map_dfr(start_indices, function(x) process_fold(x, wind_data, training_window_size, forecast_horizon))

average_lin_rmse <- mean(results$lin_rmse, na.rm = TRUE)
average_gam_rmse <- mean(results$gam_rmse, na.rm = TRUE)

print(sprintf("Average Linear Model RMSE: %f", average_lin_rmse))
print(sprintf("Average GAM Model RMSE: %f", average_gam_rmse))

fold_predictions_actuals <- data.frame(
  fold = results$fold,
  actual_spd = results$actual_spd,
  linear_prediction = results$lin_predictions,
  gam_prediction = results$gam_predictions
)

fold_predictions_actuals_long <- tidyr::pivot_longer(
  fold_predictions_actuals,
  cols = c(linear_prediction, gam_prediction),
  names_to = "Model",
  values_to = "Predicted"
) %>%
  mutate(Model = ifelse(Model == "linear_prediction", "Linear", "GAM"))

results_long <- tidyr::pivot_longer(results, cols = c(lin_rmse, gam_rmse), names_to = "Model", values_to = "RMSE")

# Calculate means for every 20 points for smoother line
smoothed_results <- results_long %>%
  group_by(Model, fold_group = ceiling(as.numeric(as.character(fold)) / 20)) %>%
  summarise(RMSE = mean(RMSE), .groups = 'drop')

# Plot with only smoothed line
ggplot(results_long, aes(x = as.numeric(as.character(fold)), y = RMSE, color = Model)) +
  geom_point(alpha = 0.3) +  # Reduced alpha for points
  geom_smooth(data = smoothed_results, aes(x = fold_group * 20 - 19, y = RMSE, group = Model), 
              se = FALSE, method = "loess", span = 0.5, size = 1) +  # Smoothing line
  scale_x_continuous(breaks = c(1, 50, 100, 150, 200, 250)) +
  labs(title = "RMSE Across Folds for Linear and GAM Models", x = "Fold", y = "RMSE") +
  theme_minimal() +
  scale_color_manual(values = c("lin_rmse" = "darkblue", "gam_rmse" = "darkorange"))

# Ensure the smoothed_results is sorted by 'fold' for proper line drawing in `geom_smooth()`
smoothed_results <- smoothed_results %>% arrange(fold_group, Model)

fold_predictions_actuals <- data.frame(
  fold = results$fold,
  actual_spd = results$actual_spd,
  linear_prediction = results$lin_predictions,
  gam_prediction = results$gam_predictions
)

fold_predictions_actuals_long <- tidyr::pivot_longer(
  fold_predictions_actuals,
  cols = c(linear_prediction, gam_prediction),
  names_to = "Model",
  values_to = "Predicted"
) %>%
  mutate(Model = ifelse(Model == "linear_prediction", "Linear", "GAM"))

ggplot(fold_predictions_actuals_long, aes(x = actual_spd, y = Predicted, color = Model)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", se = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Predicted vs. Actual Wind Speeds with LOESS", x = "Actual Wind Speed", y = "Predicted Wind Speed") +
  theme_minimal() +
  scale_color_manual(values = c("Linear" = "darkblue", "GAM" = "darkorange"))

ggplot(fold_predictions_actuals_long, aes(x = Predicted, fill = Model)) +
  geom_density(alpha = 0.7) +
  geom_density(data = fold_predictions_actuals_long, aes(x = actual_spd, y = ..density..), color = "black", fill = NA) +
  labs(title = "Density of Predicted Wind Speeds vs. Actual", x = "Wind Speed", y = "Density") +
  theme_minimal() +
  scale_fill_manual(values = c("Linear" = "darkblue", "GAM" = "darkorange"))


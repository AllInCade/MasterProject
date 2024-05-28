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
library(reshape2)
library(purrr)
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
library(DHARMa)
library(MASS)
library(TTR)
library(furrr)
library(progressr)
library(glmnet)


plan(sequential)
gc()
options(future.seed = TRUE, future.rng.onMisuse="ignore")

wind_data_ful <- read_csv("SB_Aeolus_20230519.csv")
wind_data_full <- na.omit(wind_data_ful)
rm(wind_data_ful)

wind_data_full$datetime_UTC <- ymd_hms(wind_data_full$datetime_UTC)

gc()
wind_data <- wind_data_full %>%
  filter(year(datetime_UTC) == 2021) %>%
  filter(month(datetime_UTC) >= 08) %>%
  filter(abb == 46078) %>%
  arrange(datetime_UTC) %>%
  mutate(time = row_number())

rm(wind_data_full)
gc()

# Adding a 24-hour moving average for wind speed
wind_data <- wind_data %>%
  mutate(R24 = rollmean(spd, 24, fill = NA, align = "right"))

# Adding Fourier terms for daily seasonality
wind_data <- wind_data %>%
  mutate(
    hour_of_day = hour(datetime_UTC),
    sin_hour = sin(2 * pi * hour_of_day / 24),
    cos_hour = cos(2 * pi * hour_of_day / 24)
  )


wind_data$datetime_UTC <- NULL
wind_data$lat <- NULL
wind_data$lon <- NULL
wind_data$ht <- NULL
wind_data$abb <- NULL
wind_data$name <- NULL

wind_data <- wind_data %>%
  mutate(
    spd_lag1 = lag(spd, 1),
    spd_lag2 = lag(spd, 2),
    R24 = lag(R24, 1),
    U = lag(U, 1),
    V = lag(V, 1),
    U10 = lag(U10, 1),
    V10 = lag(V10, 1),
    Utau = lag(Utau, 1),
    Vtau = lag(Vtau, 1),
    Atmp = lag(Atmp, 1),
    Wtmp = lag(Wtmp, 1),
    abar = lag(abar, 1),
    dir = lag(dir, 1),
    sin_hour = lag(sin_hour, 1),
    cos_hour = lag(cos_hour, 1)
  ) %>%
  filter(!is.na(time))  

wind_data <- na.omit(wind_data)

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

gc()

# Set parameters
total_obs <- nrow(wind_data)
initial_train_size <- floor(0.8 * total_obs)  # Use 80% for initial training
forecast_horizon <- 24

plan(multisession, workers = 6)  # Adjust for potential RAM bottleneck

# Function to update models with new data
update_models <- function(prev_gamm4, prev_gam, new_data) {
  gamm4_model <- gamm4(spd ~ s(time) + spd_lag1 + Atmp + R24 + spd_lag2,
                       data = new_data,
                       family = gaussian,
                       REML = TRUE,
                       start = prev_gamm4)
  
  gam_model <- glmmTMB(spd ~ s(time) + spd_lag1 + Wtmp + Atmp + R24 + abar,
                       disp =~ spd_lag2,
                       data = new_data,
                       family = gaussian,
                       REML = TRUE,
                       start = prev_gam)
  
  return(list(gamm4_model = gamm4_model, gam_model = gam_model))
}

# Function to process each fold
process_fold <- function(start_index, data, training_window_size, forecast_horizon, prev_models = NULL) {
  train_end <- start_index + training_window_size - 1
  test_end <- train_end + forecast_horizon
  
  train_data <- data[start_index:train_end, ]
  test_data <- data[(train_end + 1):test_end, ]
  
  if (is.null(prev_models)) {
    gamm4_model <- gamm4(spd ~ s(time) + spd_lag1 + Atmp + R24 + spd_lag2,
                         data = train_data,
                         family = gaussian,
                         REML = TRUE)
    
    gam_model <- glmmTMB(spd ~ s(time) + spd_lag1 + Wtmp + Atmp + R24 + abar,
                         disp =~ spd_lag2,
                         data = train_data,
                         family = gaussian,
                         REML = TRUE)
  } else {
    models <- update_models(prev_models$gamm4_model, prev_models$gam_model, train_data)
    gamm4_model <- models$gamm4_model
    gam_model <- models$gam_model
  }
  
  gamm4_predictions <- predict(gamm4_model$gam, newdata = test_data, type = "response")
  gam_predictions <- predict(gam_model, newdata = test_data, type = "response")
  
  gamm4_rmse <- calculate_rmse(test_data$spd, gamm4_predictions)
  gam_rmse <- calculate_rmse(test_data$spd, gam_predictions)
  
  return(list(
    gamm4_rmse = gamm4_rmse,
    gam_rmse = gam_rmse,
    actual_spd = test_data$spd,
    gamm4_predictions = gamm4_predictions,
    gam_predictions = gam_predictions,
    fold = rep(start_index, nrow(test_data)),
    gamm4_model = gamm4_model,
    gam_model = gam_model
  ))
}

# Initial training
initial_train_data <- wind_data[1:initial_train_size, ]
gamm4_initial <- gamm4(spd ~ s(time) + spd_lag1 + Atmp + R24 + spd_lag2,
                       data = initial_train_data,
                       family = gaussian,
                       REML = TRUE)

gam_initial <- glmmTMB(spd ~ s(time) + spd_lag1 + Wtmp + Atmp + R24 + abar,
                       disp =~ spd_lag2,
                       data = initial_train_data,
                       family = gaussian,
                       REML = TRUE)

initial_models <- list(gamm4_model = gamm4_initial, gam_model = gam_initial)

# Start indices for sliding window
start_indices <- seq(initial_train_size + 1, total_obs - forecast_horizon + 1, by = forecast_horizon)

# Process folds
results <- future_map_dfr(start_indices, function(x) {
  process_fold(x, wind_data, initial_train_size, forecast_horizon, initial_models)
})

# Aggregate and visualize results
error_data <- results %>%
  mutate(hour = rep(1:forecast_horizon, length(start_indices))) %>%
  pivot_longer(cols = c(gamm4_rmse, gam_rmse), names_to = "model", values_to = "rmse")

# Calculate mean RMSE for each hour
hourly_rmse <- error_data %>%
  group_by(hour, model) %>%
  summarize(mean_rmse = mean(rmse, na.rm = TRUE))

# Plot the errors
ggplot(hourly_rmse, aes(x = hour, y = mean_rmse, color = model)) +
  geom_line() +
  labs(title = "Prediction Accuracy by Hour",
       x = "Hour Ahead",
       y = "Mean RMSE",
       color = "Model") +
  theme_minimal()














































plan(sequential) # Reset back-end paralellization scheme
gc()
options(future.seed = TRUE, future.rng.onMisuse="ignore")

total_obs <- nrow(wind_data)
training_window_size <- floor(0.7*total_obs)
forecast_horizon <- 24
plan(multisession, workers = 6) # adjust workers / potential RAM bottleneck

process_fold <- function(start_index, data, training_window_size, forecast_horizon) {
  train_end <- start_index + training_window_size - 1
  test_end <- train_end + forecast_horizon
  
  train_data <- data[start_index:train_end, ]
  test_data <- data[(train_end + 1):test_end, ]
  
  gamm4_model <- gamm4(spd ~ s(time) + spd_lag1 + 
                        Atmp + R24 + spd_lag2,
                        data = train_data,
                        family = gaussian,
                        REML = TRUE)
  
  gam_model <- glmmTMB(spd ~ s(time) + spd_lag1 + 
                       Wtmp + Atmp + R24 + abar,
                       disp =~ spd_lag2,
                       data = train_data,
                       family = gaussian, 
                       REML = TRUE)
  
  gamm4_predictions <- predict(gamm4_model$gam, newdata = test_data, type = "response")
  gam_predictions <- predict(gam_model, newdata = test_data, type = "response")
  
  gamm4_rmse <- calculate_rmse(test_data$spd, gamm4_predictions)
  gam_rmse <- calculate_rmse(test_data$spd, gam_predictions)
  
  return(list(
    gamm4_rmse = gamm4_rmse,
    gam_rmse = gam_rmse,
    actual_spd = test_data$spd,
    gamm4_predictions = gamm4_predictions,
    gam_predictions = gam_predictions,
    fold = rep(start_index, nrow(test_data))  # Use 'start_index' instead of 'fold_index'
  ))
}

start_indices <- 1:(nrow(wind_data) - training_window_size - forecast_horizon + 1)
results <- future_map_dfr(start_indices, function(x) process_fold(x, wind_data, training_window_size, forecast_horizon))

average_gamm4_rmse <- mean(results$gamm4_rmse, na.rm = TRUE)
average_gam_rmse <- mean(results$gam_rmse, na.rm = TRUE)

# Create a data frame with prediction errors by hour
error_data <- results %>%
  mutate(hour = rep(1:forecast_horizon, length.out = n())) %>%
  pivot_longer(cols = c(gamm4_rmse, gam_rmse), names_to = "model", values_to = "rmse")

# Calculate mean RMSE for each hour
hourly_rmse <- error_data %>%
  group_by(hour, model) %>%
  summarize(mean_rmse = mean(rmse, na.rm = TRUE))

# Plot the errors
ggplot(hourly_rmse, aes(x = hour, y = mean_rmse, color = model)) +
  geom_line() +
  labs(title = "Prediction Accuracy by Hour",
       x = "Hour Ahead",
       y = "Mean RMSE",
       color = "Model") +
  theme_minimal()

print(sprintf("Average gamm Model RMSE: %f", average_gamm4_rmse))
print(sprintf("Average glmmTMB Model RMSE: %f", average_gam_rmse))

fold_predictions_actuals <- data.frame(
  fold = results$fold,
  actual_spd = results$actual_spd,
  gamm4_prediction = results$gamm4_predictions,
  gam_prediction = results$gam_predictions
)

fold_predictions_actuals_long <- tidyr::pivot_longer(
  fold_predictions_actuals,
  cols = c(gamm4_prediction, gam_prediction),
  names_to = "Model",
  values_to = "Predicted"
) %>%
  mutate(Model = ifelse(Model == "gamm4_prediction", "gamm4", "GAM"))

results_long <- tidyr::pivot_longer(results, cols = c(gamm4_rmse, gam_rmse), names_to = "Model", values_to = "RMSE")

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
  labs(title = "RMSE Across Folds for gamm4 and GAM Models", x = "Fold", y = "RMSE") +
  theme_minimal() +
  scale_color_manual(values = c("gamm4_rmse" = "darkblue", "gam_rmse" = "darkorange"))

# Ensure the smoothed_results is sorted by 'fold' for proper line drawing in `geom_smooth()`
smoothed_results <- smoothed_results %>% arrange(fold_group, Model)

fold_predictions_actuals <- data.frame(
  fold = results$fold,
  actual_spd = results$actual_spd,
  gamm4_prediction = results$gamm4_predictions,
  gam_prediction = results$gam_predictions
)

fold_predictions_actuals_long <- tidyr::pivot_longer(
  fold_predictions_actuals,
  cols = c(gamm4_prediction, gam_prediction),
  names_to = "Model",
  values_to = "Predicted"
) %>%
  mutate(Model = ifelse(Model == "gamm4_prediction", "gamm4", "GAM"))

ggplot(fold_predictions_actuals_long, aes(x = actual_spd, y = Predicted, color = Model)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "loess", se = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Predicted vs. Actual Wind Speeds with LOESS", x = "Actual Wind Speed", y = "Predicted Wind Speed") +
  theme_minimal() +
  scale_color_manual(values = c("gamm4" = "darkblue", "GAM" = "darkorange"))

ggplot(fold_predictions_actuals_long, aes(x = Predicted, fill = Model)) +
  geom_density(alpha = 0.7) +
  geom_density(data = fold_predictions_actuals_long, aes(x = actual_spd, y = ..density..), color = "black", fill = NA) +
  labs(title = "Density of Predicted Wind Speeds vs. Actual", x = "Wind Speed", y = "Density") +
  theme_minimal() +
  scale_fill_manual(values = c("gamm4" = "darkblue", "GAM" = "darkorange"))


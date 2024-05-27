library(quantmod)
library(mgcv)
library(glmmTMB)
library(forecast)
library(ggplot2)
library(DHARMa)
library(reshape2)
library(dplyr)
library(MASS)
library(TTR)
library(future)
library(furrr)
library(progressr)
library(glmnet)

plan(sequential)
gc()
options(future.seed = TRUE, future.rng.onMisuse="ignore")
# Load data and prepare the merged dataframe as before
japan_data <- getSymbols("^N225", from = "2020-01-01", to = "2024-04-01", auto.assign = FALSE)
euro_data <- getSymbols("^N100", from = "2020-01-01", to = "2024-04-01", auto.assign = FALSE)

# Convert to data frame and handle date indexing
japan_df <- data.frame(Date = index(japan_data), coredata(japan_data))
euro_df <- data.frame(Date = index(euro_data), coredata(euro_data))

merged_df <- merge(japan_df, euro_df, by = "Date", all = FALSE)

#merged_df$by <- NULL
merged_df <- na.locf(merged_df)
merged_df$time <- 1:nrow(merged_df)
#merged_df$date <- rownames(merged_data)

merged_df$N100_log_return <- c(NA, diff(log(merged_df$N100.Close)))
merged_df$N100_log_return_lag <- lag(merged_df$N100_log_return, 1)
merged_df$N225_log_return <- c(NA, diff(log(merged_df$N225.Close)))

# Apply lagging to other N100 columns except for those already containing "lag" in their names
for (col in grep("N100", names(merged_df), value = TRUE)) {
  if (!grepl("lag", col) && !grepl("log_return", col)) {  # Exclude log return already lagged
    lagged_col_name <- paste0(col, "_lag")
    merged_df[[lagged_col_name]] <- lag(merged_df[[col]], 1)
    merged_df[[lagged_col_name]] <- na.locf(merged_df[[lagged_col_name]], na.rm = FALSE)
  }
}

# Calculate and lag additional features
merged_df$N225_RSI <- RSI(merged_df$N225.Close, n = 20)
merged_df$N225_SMA <- SMA(merged_df$N225.Close, n = 20)
merged_df$N225_runSD <- runSD(merged_df$N225.Close, n = 20, sample = FALSE)
merged_df$N225_OC <- merged_df$N225.Open - merged_df$N225.Close

merged_df$N100_OC_lag <- lag(merged_df$N100.Close - merged_df$N100.Open, 1)
merged_df$N100_RSI_lag <- stats::lag(RSI(merged_df$N100.Close, n = 20), k = 1)
merged_df$N100_SMA_lag <- stats::lag(SMA(merged_df$N100.Close, n = 20), k = 1)
merged_df$N100_runSD_lag <- stats::lag(runSD(merged_df$N100.Close, n = 20), k = 1)

# Handle NA values more carefully to preserve lagging effect
merged_df <- na.locf(merged_df, fromLast = TRUE)  # Filling from the last to respect time order

# Decomposition and lagging seasonal and trend components
ts_data <- ts(merged_df$N100_log_return, frequency = 20)
decomposed <- stl(ts_data, s.window = "periodic", robust = TRUE)
merged_df$N100_seasonal <- stats::lag(decomposed$time.series[, "seasonal"], k = 1)
merged_df$N100_trend <- stats::lag(decomposed$time.series[, "trend"], k = 1)

# Ensure seasonal and trend components for N225 are calculated without unwanted forward fill
ts_data <- ts(merged_df$N225_log_return, frequency = 20)
decomposed <- stl(ts_data, s.window = "periodic", robust = TRUE)
merged_df$N225_seasonal <- decomposed$time.series[, "seasonal"]
merged_df$N225_trend <- decomposed$time.series[, "trend"]

# Final forward fill to handle any remaining NAs after lagging
merged_df <- na.locf(merged_df)

# Eliminate future information confounders for importance analysis
merged_df$N100.Open <- NULL
merged_df$N100.Close <- NULL
merged_df$N100.High <- NULL
merged_df$N100.Low <- NULL
merged_df$N100.Adjusted <- NULL
merged_df$N100.Volume <- NULL

rf_model <- randomForest(N100_log_return ~ ., data = merged_df, ntree = 1000, importance = TRUE)
importance_matrix <- importance(rf_model)
sorted_importance_mse <- sort(importance_matrix[, "%IncMSE"], decreasing = TRUE)
sorted_importance_purity <- sort(importance_matrix[, "IncNodePurity"], decreasing = TRUE)
print(sorted_importance_mse)
print(sorted_importance_purity)

numeric_df <- merged_df[sapply(merged_df, is.numeric)]
x <- as.matrix(numeric_df[, setdiff(names(numeric_df), "N100_log_return")])
y <- numeric_df$N100_log_return
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

total_obs <- nrow(merged_df)
training_window_size <- floor(0.8*total_obs)
forecast_horizon <- 1  
plan(multisession, workers = detectCores(logical = TRUE) - 2) # potential RAM bottleneck

process_fold <- function(start_index, data, training_window_size, forecast_horizon) {
  train_end <- start_index + training_window_size - 1
  test_end <- train_end + forecast_horizon
  
  train_data <- data[start_index:train_end, ]
  test_data <- data[(train_end + 1):test_end, ]
  
  lin_model <- lm(N100_log_return ~ time + 
                    N100_RSI_lag + N100_trend +
                    N100_log_return_lag + N225_RSI +
                    N225_log_return,
                    data = train_data)
  
  gam_model <- glmmTMB(N100_log_return ~ s(time) +
                       N100_log_return_lag + 
                       s(N100_RSI_lag) + N225_log_return +
                       s(N100_trend) + s(N225_RSI),
                       disp =~ N100.High_lag,
                       data = train_data,
                       family = gaussian(), 
                       REML = TRUE)
  
  lin_predictions <- predict(lin_model, newdata = test_data)
  gam_predictions <- predict(gam_model, newdata = test_data, type = "response")
  
  lin_rmse <- calculate_rmse(test_data$N100_log_return, lin_predictions)
  gam_rmse <- calculate_rmse(test_data$N100_log_return, gam_predictions)
  
  return(list(
    lin_rmse = lin_rmse,
    gam_rmse = gam_rmse,
    actual_log_return = test_data$N100_log_return,
    lin_predictions = lin_predictions,
    gam_predictions = gam_predictions,
    fold = rep(start_index, nrow(test_data))  # Use 'start_index' instead of 'fold_index'
  ))
}

start_indices <- 1:(nrow(merged_df) - training_window_size - forecast_horizon + 1)
results <- future_map_dfr(start_indices, function(x) process_fold(x, merged_df, training_window_size, forecast_horizon))

average_lin_rmse <- mean(results$lin_rmse, na.rm = TRUE)
average_gam_rmse <- mean(results$gam_rmse, na.rm = TRUE)

print(sprintf("Average Linear Model RMSE: %f", average_lin_rmse))
print(sprintf("Average GAM Model RMSE: %f", average_gam_rmse))

fold_predictions_actuals <- data.frame(
  fold = results$fold,
  actual_log_return = results$actual_log_return,
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
  scale_x_continuous(breaks = c(1, 50, 100, 150, 200)) +
  labs(title = "RMSE Across Folds for Linear and GAM Models", x = "Fold", y = "RMSE") +
  theme_minimal() +
  scale_color_manual(values = c("lin_rmse" = "darkblue", "gam_rmse" = "darkorange"))

# Ensure the smoothed_results is sorted by 'fold' for proper line drawing in `geom_smooth()`
smoothed_results <- smoothed_results %>% arrange(fold_group, Model)

fold_predictions_actuals <- data.frame(
  fold = results$fold,
  actual_log_return = results$actual_log_return,
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

ggplot(fold_predictions_actuals_long, aes(x = actual_log_return, y = Predicted, color = Model)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", se = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Predicted vs. Actual Log Returns with LOESS", x = "Actual Log Return", y = "Predicted Log Return") +
  theme_minimal() +
  scale_color_manual(values = c("Linear" = "darkblue", "GAM" = "darkorange"))

ggplot(fold_predictions_actuals_long, aes(x = Predicted, fill = Model)) +
  geom_density(alpha = 0.7) +
  geom_density(data = fold_predictions_actuals_long, aes(x = actual_log_return, y = ..density..), color = "black", fill = NA) +
  labs(title = "Density of Predicted Log Returns vs. Actual", x = "Log Return", y = "Density") +
  theme_minimal() +
  scale_fill_manual(values = c("Linear" = "darkblue", "GAM" = "darkorange"))

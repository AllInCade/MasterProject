library(parallel)
library(doParallel)
library(dplyr)
library(zoo)
library(lubridate)
library(readr)
library(mgcv)
library(GAM)
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

gc()
wind_data_full <- read_csv("wind_data.csv")

wind_data <- wind_data_full %>%
  filter(datetime_UTC >= ymd_hms('2007-09-10 00:00:00') & datetime_UTC <= ymd_hms('2007-10-15 23:59:59'))

wind_ts <- ts(wind_data$spd, frequency = 24) # capture diurnal pattern
sarima_model <- auto.arima(wind_ts)
residuals_data <- residuals(sarima_model)
wind_data$sarima_res <- residuals_data

wind_data <- wind_data %>%
  arrange(datetime_UTC) %>% # Make sure the data is sorted by datetime_UTC
  mutate(time = as.numeric(difftime(datetime_UTC, min(datetime_UTC), units = "hours")) + 1)

wind_data$datetime_UTC <- NULL

set.seed(123) 
wind_data <- wind_data[order(wind_data$time), ]
total_rows <- nrow(wind_data)
train_size <- floor(0.6 * total_rows)
valid_size <- floor(0.2 * total_rows)
test_size <- total_rows - train_size - valid_size 

train_data <- wind_data[1:train_size, ]
valid_data <- wind_data[(train_size + 1):(train_size + valid_size), ]
test_data <- wind_data[(train_size + valid_size + 1):total_rows, ]

calculate_rmse <- function(actual, predicted) {
  sqrt(mean((actual - predicted)^2))
}

n_folds <- 5
min_initial_obs <- floor(0.8 * train_size)
remaining_obs <- train_size - min_initial_obs
adjusted_fold_size <- floor(remaining_obs / (n_folds - 1))

Linear_rmse_values <- c()
GAM_rmse_values <- c()

pb <- progress_bar$new(
  format = "  Folding [:bar] :percent :etas",
  total = n_folds, clear = FALSE, width = 30
)

for (i in 1:n_folds) {
  
  if (i == 1) {
    current_train_end <- min_initial_obs
  } else {
    current_train_end <- min_initial_obs + (i - 1) * adjusted_fold_size
  }
  
  current_train_data <- train_data[1:current_train_end, ]
  
  if (i < n_folds) {
    current_valid_data <- train_data[(current_train_end + 1):(current_train_end + adjusted_fold_size), ]
  } else {
    # Ensure we don't exceed the bounds of the data
    final_valid_index <- min(train_size, current_train_end + 1 + adjusted_fold_size)
    current_valid_data <- train_data[(current_train_end + 1):final_valid_index, ]
  }
  
  Linear_cv <- lm(spd ~ time + U + V + Utau + Vtau + lon + lat, data = current_train_data)
  
  GAM_cv <- gam(spd ~ s(time) + s(U,V) + s(Utau,Vtau) + s(lat,lon),
                        data = current_train_data,
                        family = gaussian(link = "identity"),
                        REML = TRUE)
  
  Linear_predictions <- predict(Linear_cv, current_valid_data, type="response")
  GAM_predictions <- predict(GAM_cv, current_valid_data, type="response")
  
  Linear_rmse <- calculate_rmse(current_valid_data$spd, Linear_predictions)
  GAM_rmse <- calculate_rmse(current_valid_data$spd, GAM_predictions)
  
  Linear_rmse_values <- c(Linear_rmse_values, Linear_rmse)
  GAM_rmse_values <- c(GAM_rmse_values, GAM_rmse)
  
  cat(sprintf("Fold %d: Linear Model RMSE = %f, GAM Model RMSE = %f\n", i, Linear_rmse, GAM_rmse))
  
  pb$tick()
}

average_Linear_rmse <- mean(Linear_rmse_values)
average_GAM_rmse <- mean(GAM_rmse_values)

cat(sprintf("Average Linear Model RMSE: %f\n", average_Linear_rmse))
cat(sprintf("Average GAM Model RMSE: %f\n", average_GAM_rmse))

full_train_data <- rbind(train_data, valid_data)

start_time_Linear_train <- Sys.time()
full_model_Linear <- lm(spd ~ time + U + V + Utau + Vtau + lat + lon,
                         data = full_train_data)

end_time_Linear_train <- Sys.time()
train_duration_Linear <- end_time_Linear_train - start_time_Linear_train
start_time_Linear_pred <- Sys.time()
predictions_Linear <- predict(full_model_Linear, test_data, type = "response")
end_time_Linear_pred <- Sys.time()
pred_duration_Linear <- end_time_Linear_pred - start_time_Linear_pred
rmse_Linear <- sqrt(mean((predictions_Linear - test_data$spd)^2))

start_time_GAM_train <- Sys.time()

full_model_GAM <- gam(spd ~ s(time) + s(U,V) + s(Utau,Vtau) + s(lat,lon),
                              data = full_train_data,
                              family = gaussian(link = "identity"),
                              REML = TRUE)

end_time_GAM_train <- Sys.time()
train_duration_GAM <- end_time_GAM_train - start_time_GAM_train

start_time_GAM_pred <- Sys.time()
predictions_GAM <- predict(full_model_GAM, test_data, type = "response")
end_time_GAM_pred <- Sys.time()
pred_duration_GAM <- end_time_GAM_pred - start_time_GAM_pred
rmse_GAM <- sqrt(mean((predictions_GAM - test_data$spd)^2))

cat("Training duration for Linear: ", train_duration_Linear, "\n")
cat("Training duration for GAM: ", train_duration_GAM, "\n")

cat("Prediction duration for Linear: ", pred_duration_Linear, "\n")
cat("Prediction duration for GAM: ", pred_duration_GAM, "\n")

cat("RMSE for Linear: ", rmse_Linear, "\n")
cat("RMSE for GAM: ", rmse_GAM, "\n")

plot_data <- data.frame(
  Actual = rep(test_data$spd, 2),
  Predicted = c(predictions_Linear, predictions_GAM),
  Model = factor(rep(c("Linear", "GAM"), each = nrow(test_data)))
)

ggplot(plot_data, aes(x = Actual, y = Predicted, color = Model)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Predicted vs Actual Values", 
       x = "Actual Values", 
       y = "Predicted Values",
       color = "Model") +
  theme_minimal() +
  scale_color_manual(values = c("darkblue", "darkorange"))



ggplot(plot_data, aes(x = Predicted, fill = Model)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density of Predicted Values", 
       x = "Predicted Values",
       fill = "Model") +
  theme_minimal() +
  scale_fill_manual(values = c("Linear" = "darkorange", "GAM" = "darkblue")) +
  theme(legend.position = "right") # Include legend


actual_density <- density(full_train_data$spd)

# Convert the density object to a data frame for plotting
actual_density_df <- data.frame(Value = actual_density$x, Density = actual_density$y)

ggplot() +
  # Add the actual values density with geom_line
  geom_line(data = actual_density_df, aes(x = Value, y = Density, color = "Actual"), size = 1) +
  # Add the predicted densities
  geom_density(data = plot_data[plot_data$Model == "GAM", ], 
               aes(x = Predicted, fill = "GAM"), alpha = 0.5) +
  geom_density(data = plot_data[plot_data$Model == "Linear", ], 
               aes(x = Predicted, fill = "Linear"), alpha = 0.5) +
  scale_color_manual(values = c("Actual" = "black")) + 
  scale_fill_manual(values = c("Linear" = "darkorange", "GAM" = "darkblue")) +
  labs(title = "Density of Predicted Values vs Actual Values",
       x = "Values",
       y = "Density",
       fill = "Model",
       color = "Density Type") + 
  theme_minimal() +
  theme(legend.position = "right") +
  facet_wrap(~Model)


library(parallel)
library(doParallel)
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
gc()
wind_data_full <- read_csv("SB_Aeolus_20230519.csv")
covariates <- c("lat", "lon", "Vtau", "Utau", "U", "V", "spd", "datetime_UTC")

wind_data_full <- wind_data_full %>%
  dplyr::mutate(datetime_UTC = as.POSIXct(datetime_UTC, format="%Y-%m-%d %H:%M:%S")) %>%
  filter(datetime_UTC >= as.POSIXct("2007-08-01") & datetime_UTC <= as.POSIXct("2007-12-01")) %>%
  filter(tolower(abb) != "sm") %>%
  dplyr::mutate(numeric_time = as.numeric(difftime(datetime_UTC, as.POSIXct("1970-01-01"), units = "hours")),
         time = numeric_time - min(numeric_time)) %>%
  dplyr::mutate(across(everything(), forecast::na.interp)) %>%
  filter(spd >= 0.1) %>%
  mutate(spd_rolling_avg = rollapply(spd, width = 4, FUN = function(x) mean(x, na.rm = TRUE), 
                                     align = 'right', fill = NA, partial = TRUE)) %>%
  dplyr::select(all_of(covariates)) 

max(wind_data_full$spd)

wind_data <- wind_data_full
write_csv(wind_data, "wind_data_cleaned.csv")
rm(wind_data_full)
gc()
sample_interval <- 1  # Can be adjusted as needed
wind_data <- wind_data[seq(1, nrow(wind_data), sample_interval), ]

wind_ts <- ts(wind_data$spd, frequency = 24) # capture diurnal pattern
sarima_model <- auto.arima(wind_ts)
residuals_data <- residuals(sarima_model)
wind_data$sarima_res <- residuals_data

wind_data <- wind_data %>%
  arrange(datetime_UTC) %>% # Make sure the data is sorted by datetime_UTC
  mutate(time = as.numeric(difftime(datetime_UTC, min(datetime_UTC), units = "hours")) + 1)


set.seed(101) 
training_indices <- sample(1:nrow(wind_data), 0.8 * nrow(wind_data))
train_data <- wind_data[training_indices, ]
test_data <- wind_data[-training_indices, ]

model <- randomForest(spd ~ . - spd, data = train_data, na.action = na.omit, importance = TRUE)

predictions <- predict(model, test_data)
actuals <- test_data$spd

RMSE <- sqrt(mean((predictions - actuals)^2))
print(paste("RMSE:", RMSE))

importance_scores <- importance(model, type = 2) 
feature_names <- row.names(importance_scores)
importance_data <- data.frame(Feature = feature_names, Importance = importance_scores[, 1])
sorted_importance_data <- importance_data[order(-importance_data$Importance), ] # Sort in decreasing order
barplot(sorted_importance_data$Importance, names.arg = sorted_importance_data$Feature, 
        main = "Feature Importance in Predicting speed", las = 2, cex.names = 0.7)

linmod <- glmmTMB(spd ~ V + U + time + sarima_res, 
                  data = wind_data)

tmbmod <- glmmTMB(spd ~ s(time) + s(U) + s(V) + s(sarima_res),
                  data = wind_data, 
                  family = gaussian(link = "identity"),
                  REML = TRUE)

AIC(linmod)
AIC(tmbmod) 

set.seed(420)  
n <- nrow(wind_data)
train_size <- round(0.6 * n)
valid_size <- round(0.2 * n)
wind_data_shuffled <- wind_data[sample(nrow(wind_data)), ]
train_data <- wind_data_shuffled[1:train_size, ]
valid_data <- wind_data_shuffled[(train_size + 1):(train_size + valid_size), ]
test_data <- wind_data_shuffled[(train_size + valid_size + 1):n, ]


tmbmod2_train <- glmmTMB(spd ~ s(time) + s(U,V, k = c(5,5)) + s(Utau, Vtau, k = c(5,5)) + s(lat, lon, k = c(5,5)) + s(sarima_res),
                         data = train_data,
                         family = gaussian(link="identity"),
                         REML = TRUE)


gammod2_train <- gamm4(spd ~ s(time) + s(U,V, k = c(5,5)) + s(Utau, Vtau, k = c(5,5)) + s(lat, lon, k = c(5,5)) + s(sarima_res),
                      data = train_data, 
                      family = gaussian(link = "identity"), 
                      REML = TRUE)



AIC(tmbmod2_train)
AIC(gammod2_train$mer)


# Measure predictive power by RMSE by 5-fold CV
n_folds <- 5
folds <- cut(seq(1, nrow(valid_data)), breaks = n_folds, labels = FALSE)

# Initialize lists for RMSE
rmse_gamm <- vector("list", n_folds)
rmse_glmmTMB <- vector("list", n_folds)

# Loop through each fold
for(i in 1:n_folds) {
  # Define training and test sets for this fold
  test_indices <- which(folds == i)
  train_indices <- setdiff(seq_len(nrow(valid_data)), test_indices)
  train_fold <- valid_data[train_indices, ]
  test_fold <- valid_data[test_indices, ]
  
  model_gamm <- gamm4(spd ~ s(time) + s(U,V, k = c(5,5)) + s(Utau, Vtau, k = c(5,5)) + s(lat, lon, k = c(5,5)) + s(sarima_res),
                      data = full_train_data, 
                      family = gaussian(link = "identity"), 
                      REML = TRUE) 
  
  pred_gamm <- predict(model_gamm$gam, test_fold, type = "response")
  rmse_gamm[[i]] <- sqrt(mean((pred_gamm - test_fold$spd)^2))
  
  model_glmmTMB <- glmmTMB(spd ~ s(time) + s(U,V, k = c(5,5)) + s(Utau, Vtau, k = c(5,5)) + s(lat, lon, k = c(5,5)) + s(sarima_res),
                           disp =~ time,
                           data = full_train_data,
                           family = gaussian(link="identity"),
                           REML = TRUE)
  
  pred_glmmTMB <- predict(model_glmmTMB, test_fold, type = "response", allow.new.levels=TRUE)
  rmse_glmmTMB[[i]] <- sqrt(mean((pred_glmmTMB - test_fold$spd)^2))
}

mean_rmse_gamm <- mean(unlist(rmse_gamm), na.rm = TRUE) # Using na.rm=TRUE to handle potential NA values
mean_rmse_glmmTMB <- mean(unlist(rmse_glmmTMB), na.rm = TRUE)
cat("Average RMSE for GAMM: ", mean_rmse_gamm, "\n")
cat("Average RMSE for GLMMTMB: ", mean_rmse_glmmTMB, "\n")


full_train_data <- rbind(train_data, valid_data)

start_time_gamm_train <- Sys.time()
full_model_gamm <- gamm4(spd ~ s(time) + s(U,V, bs = "tp") + s(Utau, Vtau, bs = "tp") + s(lat, lon, bs = "tp") + s(sarima_res),
                         data = full_train_data,
                         family = gaussian(link = "identity"),
                         REML = TRUE) 

end_time_gamm_train <- Sys.time()
train_duration_gamm <- end_time_gamm_train - start_time_gamm_train
start_time_gamm_pred <- Sys.time()
predictions_gamm <- predict(full_model_gamm$gam, test_data, type = "response")
end_time_gamm_pred <- Sys.time()
pred_duration_gamm <- end_time_gamm_pred - start_time_gamm_pred
rmse_gamm <- sqrt(mean((predictions_gamm - test_data$spd)^2))

start_time_glmmTMB_train <- Sys.time()

full_model_glmmTMB <- glmmTMB(spd ~ s(time) + s(U,V, bs = "tp") + s(Utau, Vtau, bs = "tp") + s(lat, lon, bs = "tp") + s(sarima_res),
                              data = full_train_data,
                              family = gaussian(link = "identity"),
                              REML = TRUE)

end_time_glmmTMB_train <- Sys.time()
train_duration_glmmTMB <- end_time_glmmTMB_train - start_time_glmmTMB_train

start_time_glmmTMB_pred <- Sys.time()
predictions_glmmTMB <- predict(full_model_glmmTMB, test_data, type = "response")
end_time_glmmTMB_pred <- Sys.time()
pred_duration_glmmTMB <- end_time_glmmTMB_pred - start_time_glmmTMB_pred
rmse_glmmTMB <- sqrt(mean((predictions_glmmTMB - test_data$spd)^2))

cat("Training duration for GAMM: ", train_duration_gamm, "\n")
cat("Training duration for GLMMTMB: ", train_duration_glmmTMB, "\n")

cat("Prediction duration for GAMM: ", pred_duration_gamm, "\n")
cat("Prediction duration for GLMMTMB: ", pred_duration_glmmTMB, "\n")

# Print the results
cat("RMSE for GAMM: ", rmse_gamm, "\n")
cat("RMSE for GLMMTMB: ", rmse_glmmTMB, "\n")

# Correctly combine predictions for the plot
plot_data <- data.frame(
  Actual = rep(test_data$spd, 2),
  Predicted = c(predictions_gamm, predictions_glmmTMB),
  Model = factor(rep(c("GAMM", "GLMMTMB"), each = nrow(test_data)))
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
  scale_color_manual(values = c("darkorange", "darkblue"))



ggplot(plot_data, aes(x = Predicted, fill = Model)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density of Predicted Values", 
       x = "Predicted Values",
       fill = "Model") +
  theme_minimal() +
  scale_fill_manual(values = c("GAMM" = "darkorange", "GLMMTMB" = "darkblue")) +
  theme(legend.position = "right") # Include legend


actual_density <- density(full_train_data$spd)

# Convert the density object to a data frame for plotting
actual_density_df <- data.frame(Value = actual_density$x, Density = actual_density$y)

ggplot() +
  # Add the actual values density with geom_line
  geom_line(data = actual_density_df, aes(x = Value, y = Density, color = "Actual"), size = 1) +
  # Add the predicted densities
  geom_density(data = plot_data[plot_data$Model == "GLMMTMB", ], 
               aes(x = Predicted, fill = "GLMMTMB"), alpha = 0.5) +
  geom_density(data = plot_data[plot_data$Model == "GAMM", ], 
               aes(x = Predicted, fill = "GAMM"), alpha = 0.5) +
  scale_color_manual(values = c("Actual" = "black")) + 
  scale_fill_manual(values = c("GAMM" = "darkorange", "GLMMTMB" = "darkblue")) +
  labs(title = "Density of Predicted Values vs Actual Values",
       x = "Values",
       y = "Density",
       fill = "Model",
       color = "Density Type") + 
  theme_minimal() +
  theme(legend.position = "right") +
  facet_wrap(~Model)


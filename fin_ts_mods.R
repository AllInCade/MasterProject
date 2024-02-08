library(quantmod)
library(mgcv)
library(glmmTMB)
library(forecast)
library(ggplot2)
library(DHARMa)
library(reshape2)
library(dplyr)
library(MASS)

# Import time series data using quantmod
getSymbols("GME", from = "2020-12-01")

# Convert imported data to a data frame and handle missing values
GME_df <- as.data.frame(`GME`)
GME_df <- na.omit(GME_df)
# Rename columns for clarity
colnames(GME_df) <- c("open", "high", "low", "close", "volume", "adjusted")

# Create a time variable for trend analysis
GME_df$time <- 1:nrow(GME_df)

# Extract date information from the index of the xts object
GME_df$date <- index(`GME`)

# Extract year information from the date for possible seasonal analysis
GME_df$year <- format(GME_df$date, "%Y")

# Calculate log return for financial time series analysis
GME_df$log_return <- c(NA, diff(log(GME_df$close)))

# Remove NA rows created by differencing
fit_data <- na.omit(GME_df)

fit_arima <- auto.arima(fit_data$log_return, seasonal = FALSE)
fit_data$arima_residuals <- residuals(fit_arima)

set.seed(123) # For reproducibility

# Assuming 'fit_data' is your dataset
total_rows <- nrow(fit_data)
train_indices <- sample(seq_len(total_rows), size = floor(0.6 * total_rows))

train_data <- fit_data[train_indices, ]
temp_data <- fit_data[-train_indices, ]

valid_indices <- sample(seq_len(nrow(temp_data)), size = floor(0.5 * nrow(temp_data)))
valid_data <- temp_data[valid_indices, ]
test_data <- temp_data[-valid_indices, ]


# models

linmod <- glmmTMB(log_return ~ time + volume + arima_residuals, 
           data = train_data)

tmbmod<- glmmTMB(log_return ~ s(time) + s(volume) + s(arima_residuals), 
             data = train_data, REML = TRUE)



# Calculate AIC and BIC for the lm and glmmTMB models directly
aic_bic_linmod <- c(AIC = AIC(linmod), BIC = BIC(linmod))
aic_bic_tmbmod <- c(AIC = AIC(tmbmod), BIC = BIC(tmbmod))

# Combine into a data frame for nice formatting
comparison <- data.frame(
  Model = c("Linear", "GAM"),
  AIC = c(aic_bic_linmod["AIC"], aic_bic_tmbmod["AIC"]),
  BIC = c(aic_bic_linmod["BIC"], aic_bic_tmbmod["BIC"])
)

# Print the table
print(comparison)


set.seed(123)
n_folds <- 5
folds <- cut(seq(1, nrow(valid_data)), breaks = n_folds, labels = FALSE)

# Initialize lists for RMSE
rmse_lm <- vector("list", n_folds)
rmse_glmmTMB <- vector("list", n_folds)

# Loop through each fold
for(i in 1:n_folds) {
  # Define training and test sets for this fold
  test_indices <- which(folds == i)
  train_indices <- setdiff(seq_len(nrow(valid_data)), test_indices)
  train_fold <- valid_data[train_indices, ]
  test_fold <- valid_data[test_indices, ]
  
  # Train LM model
  model_lm <- glmmTMB(log_return ~ time + volume + arima_residuals, data = train_data)
  pred_lm <- predict(model_lm, test_data)
  rmse_lm[[i]] <- sqrt(mean((pred_lm - test_data$log_return)^2))
  
  # Train GLMMTMB model
  model_glmmTMB <- glmmTMB(log_return ~ s(time) + s(volume) + s(arima_residuals), data = train_data, REML = TRUE)
  pred_glmmTMB <- predict(model_glmmTMB, test_data)
  rmse_glmmTMB[[i]] <- sqrt(mean((pred_glmmTMB - test_data$log_return)^2))
}

# Calculate the average RMSE for each model
mean_rmse_lm <- mean(unlist(rmse_lm))
mean_rmse_glmmTMB <- mean(unlist(rmse_glmmTMB))

# Print the results
cat("Average RMSE for LM: ", mean_rmse_lm, "\n")
cat("Average RMSE for GLMMTMB: ", mean_rmse_glmmTMB, "\n")

full_train_data <- rbind(train_data, valid_data)

full_model_linmod <- glmmTMB(log_return ~ time + volume + arima_residuals, data = full_train_data)
full_model_tmbmod <- glmmTMB(log_return ~ s(time) + s(volume) + s(arima_residuals), data = full_train_data, REML = TRUE)

# Predict using the full model
linear_predictions <- predict(full_model_linmod, test_data)
gam_predictions <- predict(full_model_tmbmod, test_data)

plot_data <- data.frame(
  Actual = rep(test_data$log_return, 2),
  Predicted = c(predict(full_model_linmod, test_data), 
                predict(full_model_tmbmod, test_data)),
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
  scale_color_manual(values = c("darkorange", "darkblue"))


actual_density <- density(full_train_data$log_return)

actual_density_df <- data.frame(Value = actual_density$x, Density = actual_density$y)

ggplot() +
  geom_line(data = actual_density_df, aes(x = Value, y = Density, color = "Actual"), size = 1) +
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


rmse_linmod <- sqrt(mean((linear_predictions - test_data$log_return)^2))
rmse_gam <- sqrt(mean((gam_predictions - test_data$log_return)^2))
cat("Average RMSE for LM: ", rmse_linmod, "\n")
cat("Average RMSE for GLMMTMB: ", rmse_gam, "\n")
### Try different

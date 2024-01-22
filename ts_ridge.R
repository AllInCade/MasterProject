
# Load necessary libraries
library(quantmod)
library(mgcv)
library(glmmTMB)
library(forecast)
library(ggplot2)
library(DHARMa)
library(reshape2)


# Function to generate a prediction matrix for spline terms
# and re-parameterize for mixed model fitting
s2rPred <- function(sm, re, data) {
  X <- PredictMat(sm, data)   # Get prediction matrix for new data
  # Transform to random effect parameterization
  if (!is.null(re$trans.U)) X <- X %*% re$trans.U
  X <- t(t(X) * re$trans.D)
  # Re-order columns according to random effect re-ordering
  X[, re$rind] <- X[, re$pen.ind != 0] 
  # Re-order penalization index in the same way  
  pen.ind <- re$pen.ind; pen.ind[re$rind] <- pen.ind[pen.ind > 0]
  # Start return object
  r <- list(rand = list(), Xf = X[, which(re$pen.ind == 0), drop = FALSE])
  for (i in 1:length(re$rand)) { # Loop over random effect matrices
    r$rand[[i]] <- X[, which(pen.ind == i), drop = FALSE]
    attr(r$rand[[i]], "s.label") <- attr(re$rand[[i]], "s.label")
  }
  names(r$rand) <- names(re$rand)
  r
}

# Import time series data using quantmod
getSymbols("NIO", from = "2020-01-01")

# Convert imported data to a data frame and handle missing values
NIO_df <- as.data.frame(`NIO`)
NIO_df <- na.omit(NIO_df)
# Rename columns for clarity
colnames(NIO_df) <- c("open", "high", "low", "log_return", "volume", "adjusted")

# Create a time variable for trend analysis
NIO_df$time <- 1:nrow(NIO_df)

# Extract date information from the index of the xts object
NIO_df$date <- index(`NIO`)

# Extract year information from the date for possible seasonal analysis
NIO_df$year <- format(NIO_df$date, "%Y")

# Calculate log return for financial time series analysis
NIO_df$log_return <- c(NA, diff(log(NIO_df$log_return)))

# Remove NA rows created by differencing
fit_data <- na.omit(NIO_df)

calculate_n_knots <- function(smooth_term) {
  smooth_term[["bs.dim"]] - smooth_term[["null.space.dim"]]
}

# Generate smooth terms for 'time' using mgcv's smoothCon function
sm_time <- mgcv::smoothCon(s(time, k = 10, bs="cs"), absorb.cons = TRUE, data = fit_data)[[1]]
re_time <- mgcv::smooth2random(sm_time, "", type = 2)
pred_matrix_time <- s2rPred(sm_time, re_time, data = fit_data)
fit_data$Xf_time <- pred_matrix_time$Xf
fit_data$Xr_time <- pred_matrix_time$rand[[1]]

# Generate smooth terms for 'volume' in a similar fashion
sm_volume <- mgcv::smoothCon(s(volume, k = 10, bs="cs"), absorb.cons = TRUE, data = fit_data)[[1]]
re_volume <- mgcv::smooth2random(sm_volume, "", type = 2)
pred_matrix_volume <- s2rPred(sm_volume, re_volume, data = fit_data)
fit_data$Xf_volume <- pred_matrix_volume$Xf
fit_data$Xr_volume <- pred_matrix_volume$rand[[1]]


n_knots_time <- calculate_n_knots(sm_time)
n_knots_volume <- calculate_n_knots(sm_volume)


# Fit an ARIMA model to the log returns and extract residuals
# These residuals will be used in the spline regression model
fit_arima <- auto.arima(fit_data$log_return, seasonal = FALSE)
fit_data$arima_residuals <- residuals(fit_arima)




###############################################
########    Set Up Ridge Penalty      #########
###############################################



augment_data <- function(data, lambda, n_knots_list) {
  augmented_matrices <- lapply(names(n_knots_list), function(name) {
    augment_design_matrix(data[[paste0("Xr_", name)]], lambda)
  })
  
  y_augmented <- augment_response(data$degrees_from_mean, max(sapply(augmented_matrices, nrow)))
  
  augmented_data <- data.frame(y_augmented)
  for (name in names(n_knots_list)) {
    augmented_data[[paste0("Xr_", name)]] <- I(matrix(unlist(augmented_matrices[[name]]), 
                                                      ncol = n_knots_list[[name]], byrow = TRUE))
  }
  augmented_data
}




# Function to perform time-based cross-validation for lambda selection
cross_validate_lambda <- function(fit_data, lambda_values, k_folds) {
  cv_errors <- numeric(length(lambda_values))
  n <- nrow(fit_data)
  fold_size <- floor(n / k_folds)
  
  for (i in seq_along(lambda_values)) {
    lambda <- lambda_values[i]
    rmse_values <- numeric(k_folds)
    
    for (fold in 1:k_folds) {
      # Calculate indices for the validation set in a time-ordered manner
      val_start <- (fold - 1) * fold_size + 1
      val_end <- fold * fold_size
      val_end <- ifelse(fold == k_folds, n, val_end)  # Adjust for last fold
      
      train_set_augmented <- augment_data(train_set, lambda, 
                                          list(year = n_knots_year, month = n_knots_month, 
                                               max_temp = n_knots_max_temp, min_temp = n_knots_min_temp, 
                                               lat_lon = n_knots_lat_lon))
      validation_set_augmented <- augment_data(validation_set, lambda, 
                                               list(year = n_knots_year, month = n_knots_month, 
                                                    max_temp = n_knots_max_temp, min_temp = n_knots_min_temp, 
                                                    lat_lon = n_knots_lat_lon))
      
      # Modify the augment_data function call in the cross_validate_lambda function
      train_set_augmented <- augment_data(train_set, lambda, n_knots_time, n_knots_volume)
      validation_set_augmented <- augment_data(validation_set, lambda, n_knots_time, n_knots_volume)
      
      # Fit the model on the augmented training set
      ccv_model <- glmmTMB(formula = y_augmented ~ Xr_year + Xr_month + Xr_max_temp + Xr_min_temp + Xr_lat_lon, 
                           data = train_set_augmented, 
                           family = gaussian(link = "identity"))
      
      
      # Evaluate on the augmented validation set and store the RMSE
      preds <- predict(cv_model, newdata = validation_set_augmented, re.form = NA)  # get predictions
      actual <- validation_set_augmented$log_return  # actual values
      rmse_values[fold] <- sqrt(mean((preds - actual) ^ 2))  # calculate RMSE
      
    }
    
    cv_errors[i] <- mean(rmse_values)
  }
  lambda_values[which.min(cv_errors)]
}



# Define a range of lambda values for cross-validation
lambda_values <- seq(0.1, 10, by = 0.5)

# Perform cross-validation to select lambda
optimal_lambda <- cross_validate_lambda(fit_data, lambda_values, k_folds = 5)

# Augment the data with the selected lambda
fit_data_augmented <- augment_data(fit_data, optimal_lambda, n_knots_time, n_knots_volume)




#####################################################
######    Data Splitting and Retraining       #######
#####################################################



# Root Mean Squared Error function
rmse <- function(y_actual, y_pred) {
  sqrt(mean((y_actual - y_pred)^2))
}

set.seed(10)
# Set the proportion of data to use for training
train_proportion <- 0.7

# Calculate the number of observations to use for training
train_size <- floor(train_proportion * nrow(fit_data))

# Create training and testing datasets while preserving the time series order
train_data <- fit_data[1:train_size, ]
test_data <- fit_data[(train_size + 1):nrow(fit_data), ]

# Calculate the number of knots for time and volume for augmentation
n_knots_time <- calculate_n_knots(sm_time)  # Assuming calculate_n_knots is defined as before
n_knots_volume <- calculate_n_knots(sm_volume)

# Augment the training data with optimal_lambda
train_data_augmented <- augment_data(train_data, optimal_lambda, n_knots_time, n_knots_volume)

# Start timing for tmb1
start_time <- Sys.time()
tmb1_train <- glmmTMB(formula = log_return ~ Xr_time + Xr_volume + arima_residuals,
                      data = train_data_augmented, 
                      family = gaussian(link="identity"))
end_time <- Sys.time()
time_tmb1 <- end_time - start_time
cat("Time taken for tmb1: ", time_tmb1, "seconds\n")

# Start timing for tmb2
start_time <- Sys.time()
tmb2_train <- glmmTMB(formula = log_return ~ Xr_time + Xr_volume + arima_residuals,
                      data = train_data, 
                      family = gaussian(link="identity"))
end_time <- Sys.time()
time_tmb2 <- end_time - start_time
cat("Time taken for tmb2: ", time_tmb2, "seconds\n")

# Start timing for tmb3
start_time <- Sys.time()
tmb3_train <- glmmTMB(log_return ~ s(time, k = 10) + s(volume, k = 10) + arima_residuals,
                      data = train_data, 
                      family = gaussian(link="identity"))
end_time <- Sys.time()
time_tmb3 <- end_time - start_time
cat("Time taken for tmb3: ", time_tmb3, "seconds\n")




#################################################
######       Cross-Validation          ##########
#################################################



# Function to calculate Mean Absolute Error
calculate_mae <- function(model, data, response_col) {
  preds <- predict(model, newdata = data, re.form = NA)
  return(mean(abs(data[[response_col]] - preds)))
}

# Function to calculate R-squared
calculate_r_squared <- function(model, data, response_col) {
  preds <- predict(model, newdata = data, re.form = NA)
  return(cor(data[[response_col]], preds)^2)
}

# Function to calculate Mean Absolute Percentage Error
calculate_mape <- function(model, data, response_col) {
  preds <- predict(model, newdata = data, re.form = NA)
  actuals <- data[[response_col]]
  nonzero_indices <- which(actuals != 0)  # Indices where actual values are not zero
  return(mean(abs((actuals[nonzero_indices] - preds[nonzero_indices]) / actuals[nonzero_indices]), na.rm = TRUE))
}


# Wrapper function for model fitting
fit_model <- function(formula, data, family, augment = FALSE, lambda = NULL, n_knots_time = NULL, n_knots_volume = NULL) {
  if (augment && !is.null(lambda)) {
    data <- augment_data(data, lambda, n_knots_time, n_knots_volume)
  }
  return(glmmTMB(formula = formula, data = data, family = family))
}

# Wrapper function for calculating RMSE
calculate_rmse <- function(model, data, response_col) {
  preds <- predict(model, newdata = data, re.form = NA)
  return(rmse(data[[response_col]], preds))
}

# Rolling window cross-validation
set.seed(10)
k <- 10
n <- nrow(train_data)
window_size <- floor(n / k)  # Adjust window size as needed

# Initialize vectors to hold RMSEs and timing
cv_rmse <- matrix(nrow = k, ncol = 3)
timing <- matrix(nrow = k, ncol = 3)
cv_mae <- matrix(nrow = k, ncol = 3)
cv_r_squared <- matrix(nrow = k, ncol = 3)
cv_mape <- matrix(nrow = k, ncol = 3)


# List to store dataframes from each fold
all_folds_data <- list()  

# Loop for cross-validation
for(i in 1:k) {
  train_end <- i * window_size
  if (i == k) {
    train_end <- n  # Ensure all data is used in the last fold
    valid_set <- train_data[(train_end - window_size + 1):n, ]  # Adjust validation set for the last fold
  } else {
    valid_set <- train_data[(train_end + 1):min(n, train_end + window_size), ]
  }
  
  train_set <- train_data[1:train_end, ]
  
  # Fit models
  tryCatch({
    # Fit tmb1 on augmented training data
    start_time <- Sys.time()
    tmb1_model <- fit_model(log_return ~ Xr_time + Xr_volume + arima_residuals, train_set, gaussian(link="identity"), TRUE, optimal_lambda, n_knots_time, n_knots_volume)
    timing[i, 1] <- Sys.time() - start_time
    # Augment valid_set before making predictions for tmb1
    valid_set_augmented <- augment_data(valid_set, optimal_lambda, n_knots_time, n_knots_volume)
    cv_rmse[i, 1] <- calculate_rmse(tmb1_model, valid_set_augmented, "log_return")
    # Calculate additional metrics
    cv_mae[i, 1] <- calculate_mae(tmb1_model, valid_set_augmented, "log_return")
    cv_r_squared[i, 1] <- calculate_r_squared(tmb1_model, valid_set_augmented, "log_return")
    cv_mape[i, 1] <- calculate_mape(tmb1_model, valid_set_augmented, "log_return")
    
  }, error = function(e) {
    print(paste("Error in fitting tmb1 model in fold", i, ":", e$message))
    timing[i, 1] <- NA
    cv_rmse[i, 1] <- NA
  })
  # Model tmb2
  tryCatch({
    start_time <- Sys.time()
    tmb2_model <- fit_model(log_return ~ Xr_time + Xr_volume + arima_residuals, train_set, gaussian(link="identity"))
    timing[i, 2] <- Sys.time() - start_time
    cv_rmse[i, 2] <- calculate_rmse(tmb2_model, valid_set, "log_return")
    cv_mae[i, 2] <- calculate_mae(tmb2_model, valid_set, "log_return")
    cv_r_squared[i, 2] <- calculate_r_squared(tmb2_model, valid_set, "log_return")
    cv_mape[i, 2] <- calculate_mape(tmb2_model, valid_set, "log_return")
  }, error = function(e) {
    print(paste("Error in fitting tmb2 model in fold", i, ":", e$message))
    timing[i, 2] <- NA
    cv_rmse[i, 2] <- NA
  })
  
  # Model tmb3
  tryCatch({
    start_time <- Sys.time()
    tmb3_model <- fit_model(log_return ~ s(time, k = 10) + s(volume, k = 10) + arima_residuals, train_set, gaussian(link="identity"))
    timing[i, 3] <- Sys.time() - start_time
    cv_rmse[i, 3] <- calculate_rmse(tmb3_model, valid_set, "log_return")
    cv_mae[i, 3] <- calculate_mae(tmb3_model, valid_set, "log_return")
    cv_r_squared[i, 3] <- calculate_r_squared(tmb3_model, valid_set, "log_return")
    cv_mape[i, 3] <- calculate_mape(tmb3_model, valid_set, "log_return")
  }, error = function(e) {
    print(paste("Error in fitting tmb3 model in fold", i, ":", e$message))
    timing[i, 3] <- NA
    cv_rmse[i, 3] <- NA
  })
  

  # Calculate residuals and predictions
  residuals_tmb1 <- valid_set$log_return - predict(tmb1_model, newdata = valid_set_augmented, re.form = NA)
  residuals_tmb2 <- valid_set$log_return - predict(tmb2_model, newdata = valid_set, re.form = NA)
  residuals_tmb3 <- valid_set$log_return - predict(tmb3_model, newdata = valid_set, re.form = NA)
  predictions_tmb1 <- predict(tmb1_model, newdata = valid_set_augmented, re.form = NA)
  predictions_tmb2 <- predict(tmb2_model, newdata = valid_set, re.form = NA)
  predictions_tmb3 <- predict(tmb3_model, newdata = valid_set, re.form = NA)
  
  # Creating a dataframe for each fold
  fold_data_tmb1 <- data.frame(Model = "tmb1", Predicted = predictions_tmb1, Residuals = residuals_tmb1)
  fold_data_tmb2 <- data.frame(Model = "tmb2", Predicted = predictions_tmb2, Residuals = residuals_tmb2)
  fold_data_tmb3 <- data.frame(Model = "tmb3", Predicted = predictions_tmb3, Residuals = residuals_tmb3)
  
  # Combine fold data into the list
  all_folds_data[[i]] <- rbind(fold_data_tmb1, fold_data_tmb2, fold_data_tmb3)
}

# Combining all folds data into a single dataframe
residuals_data <- do.call(rbind, all_folds_data)


# Calculate average RMSE and total timing
avg_rmse <- colMeans(cv_rmse)
total_timing <- colSums(timing)
avg_mae <- colMeans(cv_mae)
avg_r_squared <- colMeans(cv_r_squared)
avg_mape <- colMeans(cv_mape)



##################################################
#######        Learning curves          ##########
##################################################



# Initialize
max_train_size <- nrow(train_data)
train_sizes <- seq(200, max_train_size, by = 50)

# Initialize arrays for RMSE storage
rmse_train_tmb1 <- numeric(length(train_sizes))
rmse_train_tmb2 <- numeric(length(train_sizes))
rmse_train_tmb3 <- numeric(length(train_sizes))

rmse_valid_tmb1 <- numeric(length(train_sizes))
rmse_valid_tmb2 <- numeric(length(train_sizes))
rmse_valid_tmb3 <- numeric(length(train_sizes))

# Initialize vectors to store timing for each model
timing_tmb1 <- numeric(length(train_sizes))
timing_tmb2 <- numeric(length(train_sizes))
timing_tmb3 <- numeric(length(train_sizes))

# Main loop for training size subsets
for (i in seq_along(train_sizes)) {
  size <- train_sizes[i]
  
  # Use contiguous blocks of data up to 'size' for the training data
  subset_train_data <- train_data[1:size, ]
  # The validation set is the data immediately following the training data
  valid_set <- train_data[(size + 1):max_train_size, ]
  
  # Augment the subset for the tmb1 model
  subset_train_data_augmented <- augment_data(subset_train_data, optimal_lambda, n_knots_time, n_knots_volume)
  
  # Fit the models on the training subset and record timing
  start_time <- Sys.time()
  tmb1_train <- glmmTMB(log_return ~ Xr_time + Xr_volume + arima_residuals, 
                        data = subset_train_data_augmented, 
                        family = gaussian(link="identity"))
  end_time <- Sys.time()
  timing_tmb1[i] <- as.numeric(end_time - start_time, units = "secs")
  
  # Repeat fitting for tmb2 and tmb3 models
  start_time <- Sys.time()
  tmb2_train <- glmmTMB(log_return ~ Xr_time + Xr_volume + arima_residuals, 
                        data = subset_train_data, 
                        family = gaussian(link="identity"))
  end_time <- Sys.time()
  timing_tmb2[i] <- as.numeric(end_time - start_time, units = "secs")
  
  start_time <- Sys.time()
  tmb3_train <- glmmTMB(log_return ~ s(time, k = 10) + s(volume, k = 10) + arima_residuals, 
                        data = subset_train_data, 
                        family = gaussian(link="identity"))
  end_time <- Sys.time()
  timing_tmb3[i] <- as.numeric(end_time - start_time, units = "secs")
  
  
  # Augment valid_set before making predictions for tmb1
  valid_set_augmented <- augment_data(valid_set, optimal_lambda, n_knots_time, n_knots_volume)
  
  # Make predictions on the training and validation set
  pred_tmb1_train <- predict(tmb1_train, newdata = subset_train_data_augmented)
  pred_tmb1_valid <- predict(tmb1_train, newdata = valid_set_augmented)
  
  pred_tmb2_train <- predict(tmb2_train, newdata = subset_train_data)
  pred_tmb2_valid <- predict(tmb2_train, newdata = valid_set)
  
  pred_tmb3_train <- predict(tmb3_train, newdata = subset_train_data)
  pred_tmb3_valid <- predict(tmb3_train, newdata = valid_set)
  
  # Calculate RMSE on the training and validation set
  rmse_train_tmb1[i] <- rmse(subset_train_data$log_return, pred_tmb1_train)
  rmse_valid_tmb1[i] <- rmse(valid_set$log_return, pred_tmb1_valid)
  
  rmse_train_tmb2[i] <- rmse(subset_train_data$log_return, pred_tmb2_train)
  rmse_valid_tmb2[i] <- rmse(valid_set$log_return, pred_tmb2_valid)
  
  rmse_train_tmb3[i] <- rmse(subset_train_data$log_return, pred_tmb3_train)
  rmse_valid_tmb3[i] <- rmse(valid_set$log_return, pred_tmb3_valid)
}


# Calculate total time for each model after the loop
total_time_tmb1 <- sum(timing_tmb1)
total_time_tmb2 <- sum(timing_tmb2)
total_time_tmb3 <- sum(timing_tmb3)

# Print the total timing results
cat("Total time for tmb1: ", total_time_tmb1, "seconds\n")
cat("Total time for tmb2: ", total_time_tmb2, "seconds\n")
cat("Total time for tmb3: ", total_time_tmb3, "seconds\n")


for (i in seq_along(train_sizes)) {
  cat("Training Size:", train_sizes[i], "\n")
  cat("RMSE on Training Data - tmb1:", rmse_train_tmb1[i], ", tmb2:", rmse_train_tmb2[i], ", tmb3:", rmse_train_tmb3[i], "\n")
  cat("RMSE on Validation Data - tmb1:", rmse_valid_tmb1[i], ", tmb2:", rmse_valid_tmb2[i], ", tmb3:", rmse_valid_tmb3[i], "\n")
  cat("\n")
}



# Prepare data for RMSE plots
rmse_data_training <- data.frame(
  TrainingSize = rep(train_sizes, 3),
  RMSE = c(rmse_train_tmb1, rmse_train_tmb2, rmse_train_tmb3),
  Model = rep(c("tmb1", "tmb2", "tmb3"), each = length(train_sizes))
)

rmse_data_validation <- data.frame(
  TrainingSize = rep(train_sizes, 3),
  RMSE = c(rmse_valid_tmb1, rmse_valid_tmb2, rmse_valid_tmb3),
  Model = rep(c("tmb1", "tmb2", "tmb3"), each = length(train_sizes))
)

# Prepare data for timing plot
timing_data <- data.frame(
  TrainingSize = rep(train_sizes, 3),
  Time = c(timing_tmb1, timing_tmb2, timing_tmb3),
  Model = rep(c("tmb1", "tmb2", "tmb3"), each = length(train_sizes))
)



ggplot(timing_data, aes(x = TrainingSize, y = Time, color = Model)) +
  geom_line(size = 1) +
  geom_point() +
  labs(title = "Model Training Time", x = "Training Size", y = "Time (seconds)", color = "Model") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Combine training and validation data for a single plot
combined_rmse_data <- rbind(
  transform(rmse_data_training, DataSet = "Training"),
  transform(rmse_data_validation, DataSet = "Validation")
)

ggplot(combined_rmse_data, aes(x = TrainingSize, y = RMSE, color = Model, linetype = DataSet)) +
  geom_line() +
  geom_smooth(method = "loess", se = FALSE) +
  labs(title = "Learning Curves (Training and Validation)", x = "Training Size", y = "RMSE", color = "Model", linetype = "Data Set") +
  theme_minimal() +
  theme(legend.position = "bottom")


ggplot(rmse_data_training, aes(x = TrainingSize, y = RMSE, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~Model, scales = "free_y") +
  labs(title = "Training RMSE by Model", x = "Training Size", y = "RMSE") +
  theme_minimal() +
  theme(legend.position = "none")


ggplot(rmse_data_validation, aes(x = TrainingSize, y = RMSE, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~Model, scales = "free_y") +
  labs(title = "Validation RMSE by Model", x = "Training Size", y = "RMSE") +
  theme_minimal() +
  theme(legend.position = "none")





######################################################
###  Statistics and Evaluations on Training Data   ###
######################################################



summary(tmb1_train)
summary(tmb2_train)
summary(tmb3_train)


aic_tmb1 <- AIC(tmb1_train)
bic_tmb1 <- BIC(tmb1_train)

aic_tmb2 <- AIC(tmb2_train)
bic_tmb2 <- BIC(tmb2_train)

aic_tmb3 <- AIC(tmb3_train)
bic_tmb3 <- BIC(tmb3_train)

print(paste('AIC for TMB Model 1: ', aic_tmb1))
print(paste('BIC for TMB Model 1: ', bic_tmb1))

print(paste('AIC for TMB Model 2: ', aic_tmb2))
print(paste('BIC for TMB Model 2: ', bic_tmb2))


print(paste('AIC for glmmTMB Model 1: ', aic_tmb3))
print(paste('BIC for glmmTMB Model 1: ', bic_tmb3))


# Residuals
residuals_tmb1 <- residuals(tmb1_train)
residuals_tmb2 <- residuals(tmb2_train)
residuals_tmb3 <- residuals(tmb3_train)

# Residual Plots
plot(residuals_tmb1, main = "Residuals for tmb1")
plot(residuals_tmb2, main = "Residuals for tmb2")
plot(residuals_tmb3, main = "Residuals for tmb3")

# Autocorrelation check 
acf(residuals(tmb1_train), main = "ACF of Residuals for tmb1")
acf(residuals(tmb2_train), main = "ACF of Residuals for tmb2")
acf(residuals(tmb3_train), main = "ACF of Residuals for tmb3")


# Simulate residuals and create DHARMa residual plots
simulation_tmb1 <- simulateResiduals(fittedModel = tmb1_train)
plot(simulation_tmb1)

simulation_tmb2 <- simulateResiduals(fittedModel = tmb2_train)
plot(simulation_tmb2)

simulation_tmb3 <- simulateResiduals(fittedModel = tmb3_train)
plot(simulation_tmb3)

# Perform and print dispersion tests using the already simulated residuals
dispersion_test_tmb1 <- testDispersion(simulation_tmb1)
print(dispersion_test_tmb1)

dispersion_test_tmb2 <- testDispersion(simulation_tmb2)
print(dispersion_test_tmb2)

dispersion_test_tmb3 <- testDispersion(simulation_tmb3)
print(dispersion_test_tmb3)



# Make in-sample predictions and compute RMSE
pred_tmb1_train <- predict(tmb1_train, newdata = train_data_augmented)
rmse_tmb1_train <- rmse(train_data_augmented$log_return, pred_tmb1_train)

pred_tmb2_train <- predict(tmb2_train, newdata = train_data)
rmse_tmb2_train <- rmse(train_data$log_return, pred_tmb2_train)

pred_tmb3_train <- predict(tmb3_train, newdata = train_data)
rmse_tmb3_train <- rmse(train_data$log_return, pred_tmb3_train)

# Display the RMSE for each model
print(paste("In-sample RMSE for tmb1_train: ", rmse_tmb1_train))
print(paste("In-sample RMSE for tmb2_train: ", rmse_tmb2_train))
print(paste("In-sample RMSE for tmb3_train: ", rmse_tmb3_train))




########################################################
###  Statistics and Evaluations on Validation Data   ###
########################################################



# Print results
cat("Average RMSE for tmb1 across", k, "folds: ", avg_rmse[1], "\n")
cat("Average RMSE for tmb2 across", k, "folds: ", avg_rmse[2], "\n")
cat("Average RMSE for tmb3 across", k, "folds: ", avg_rmse[3], "\n")
cat("Total time for tmb1 across", k, "folds: ", total_timing[1], "seconds\n")
cat("Total time for tmb2 across", k, "folds: ", total_timing[2], "seconds\n")
cat("Total time for tmb3 across", k, "folds: ", total_timing[3], "seconds\n")
cat("Average MAE for tmb1 across", k, "folds", avg_mae[1], "\n")
cat("Average MAE for tmb2 across", k, "folds", avg_mae[2], "\n")
cat("Average MAE for tmb3 across", k, "folds", avg_mae[3], "\n")
cat("Average r-squared for tmb1 across", k, "folds", avg_r_squared[1], "\n")
cat("Average r-squared for tmb2 across", k, "folds", avg_r_squared[2], "\n")
cat("Average r-squared for tmb3 across", k, "folds", avg_r_squared[3], "\n")
cat("Average mape for tmb1 across", k, "folds", avg_mape[1], "\n")
cat("Average mape for tmb2 across", k, "folds", avg_mape[2], "\n")
cat("Average mape for tmb3 across", k, "folds", avg_mape[3], "\n")



#######################################
####       Visuals and Plots       ####
#######################################



# Preparing data for RMSE plotting
rmse_data_cv <- data.frame(
  Fold = rep(1:k, 3),
  RMSE = c(cv_rmse[,1], cv_rmse[,2], cv_rmse[,3]),
  Model = rep(c("tmb1", "tmb2", "tmb3"), each = k)
)

ggplot(rmse_data_cv, aes(x = Fold, y = RMSE, color = Model)) +
  geom_line(size = 1.2, alpha = 0.8) +
  geom_point(size = 3, shape = 21, fill = "white") +
  scale_color_brewer(palette = "Set1") +
  labs(title = "RMSE Across Cross-Validation Folds", x = "Fold", y = "RMSE") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 12),
    panel.grid.major = element_line(color = "gray80"),
    panel.grid.minor = element_blank()
  )



# Preparing data for MAE plotting
mae_data_cv <- data.frame(
  Fold = rep(1:k, 3),
  MAE = c(cv_mae[,1], cv_mae[,2], cv_mae[,3]),
  Model = rep(c("tmb1", "tmb2", "tmb3"), each = k)
)

ggplot(mae_data_cv, aes(x = Fold, y = MAE, color = Model)) +
  geom_line(size = 1.2, alpha = 0.8) +
  geom_point(size = 3, shape = 21, fill = "white") +
  scale_color_brewer(palette = "Set1") +
  labs(title = "MAE Across Cross-Validation Folds", x = "Fold", y = "MAE") +
  theme_classic() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 12),
    panel.grid.major = element_line(color = "gray80")
  )

# Preparing data for R-squared plotting
rsq_data_cv <- data.frame(
  Fold = rep(1:k, 3),
  R_squared = c(cv_r_squared[,1], cv_r_squared[,2], cv_r_squared[,3]),
  Model = rep(c("tmb1", "tmb2", "tmb3"), each = k)
)

ggplot(rsq_data_cv, aes(x = Fold, y = R_squared, color = Model)) +
  geom_line(size = 1.2, alpha = 0.8) +
  geom_point(size = 3, shape = 21, fill = "white") +
  scale_color_brewer(palette = "Set1") +
  labs(title = "R-squared Across Cross-Validation Folds", x = "Fold", y = "R-squared") +
  theme_light() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 12),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  )



# Preparing data for R-squared plotting
mape_data_cv <- data.frame(
  Fold = rep(1:k, 3),
  MAPE = c(cv_mape[,1], cv_mape[,2], cv_mape[,3]),
  Model = rep(c("tmb1", "tmb2", "tmb3"), each = k)
)

ggplot(mape_data_cv, aes(x = Fold, y = MAPE, color = Model)) +
  geom_line(size = 1.2, alpha = 0.8) +
  geom_point(size = 3, shape = 21, fill = "white") +
  scale_color_brewer(palette = "Set1") +
  labs(title = "MAPE Across Cross-Validation Folds", x = "Fold", y = "MAPE") +
  theme_light() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 12),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  )



# Preparing data for timing plotting
timing_data_cv <- data.frame(
  Model = c("tmb1", "tmb2", "tmb3"),
  TotalTime = total_timing
)

ggplot(timing_data_cv, aes(x = Model, y = TotalTime, fill = Model)) +
  geom_bar(stat = "identity") +
  labs(title = "Total Model Training Time Across All Folds", x = "Model", y = "Total Time (seconds)") +
  theme_minimal() +
  theme(legend.position = "none")



ggplot(residuals_data, aes(x = Predicted, y = Residuals, color = Model)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Residuals Plot", x = "Predicted Values", y = "Residuals") +
  facet_wrap(~Model) +
  theme_minimal()



####################################
###   Fit Models on Full Data    ###
####################################


# Assuming train_data and valid_data are combined to form the full dataset
full_data <- rbind(train_set, valid_set)

# Fit tmb1 on augmented full dataset
full_data_augmented <- augment_data(full_data, optimal_lambda, n_knots_time, n_knots_volume)
tmb1_full <- glmmTMB(log_return ~ Xr_time + Xr_volume + arima_residuals, 
                     data = full_data_augmented, 
                     family = gaussian(link="identity"))

# Fit tmb2 on full dataset
tmb2_full <- glmmTMB(log_return ~ Xr_time + Xr_volume + arima_residuals, 
                     data = full_data, 
                     family = gaussian(link="identity"))

# Fit tmb3 on full dataset
tmb3_full <- glmmTMB(log_return ~ s(time, k = 10) + s(volume, k = 10) + arima_residuals, 
                     data = full_data, 
                     family = gaussian(link="identity"))


# Calculate RMSE 

# Make predictions on the test set
pred_tmb1_test <- predict(tmb1_full, newdata = test_data, re.form = NA)
pred_tmb2_test <- predict(tmb2_full, newdata = test_data, re.form = NA)
pred_tmb3_test <- predict(tmb3_full, newdata = test_data, re.form = NA)

# Calculate RMSE for each model on test data
rmse_tmb1_test <- rmse(test_data$log_return, pred_tmb1_test)
rmse_tmb2_test <- rmse(test_data$log_return, pred_tmb2_test)
rmse_tmb3_test <- rmse(test_data$log_return, pred_tmb3_test)

# Print RMSE results
cat("RMSE for tmb1 on test data: ", rmse_tmb1_test, "\n")
cat("RMSE for tmb2 on test data: ", rmse_tmb2_test, "\n")
cat("RMSE for tmb3 on test data: ", rmse_tmb3_test, "\n")




#########################################
###    Plots and Visuals Test Data    ###
#########################################


# Preparing data for plotting
test_data_predictions <- data.frame(
  Actual = rep(test_data$log_return, 3),
  Predicted = c(pred_tmb1_test, pred_tmb2_test, pred_tmb3_test),
  Model = factor(rep(c("tmb1", "tmb2", "tmb3"), each = nrow(test_data)))
)

ggplot(test_data_predictions, aes(x = Actual, y = Predicted, color = Model)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  geom_smooth(method = "loess", se = TRUE, aes(fill = Model), alpha = 0.3) +  # Add loess smooth line with CI
  labs(title = "Predicted vs Actual Values on Test Data", x = "Actual", y = "Predicted") +
  theme_minimal() +
  theme(legend.position = "bottom")



# Calculating residuals
test_data_predictions$Residuals = test_data_predictions$Actual - test_data_predictions$Predicted

ggplot(test_data_predictions, aes(x = Predicted, y = Residuals, color = Model)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(title = "Residuals on Test Data", x = "Predicted", y = "Residuals") +
  theme_minimal() +
  theme(legend.position = "bottom")



ggplot(test_data_predictions, aes(x = Predicted, fill = Model)) +
  geom_density(alpha = 0.7) +
  labs(title = "Density of Predictions on Test Data", x = "Predicted", y = "Density") +
  theme_minimal() +
  theme(legend.position = "top") +
  facet_wrap(~Model)



ggplot(test_data_predictions, aes(x = Residuals, fill = Model)) +
  geom_histogram(bins = 30, alpha = 0.7) +
  labs(title = "Histogram of Prediction Errors", x = "Error", y = "Count") +
  theme_minimal() +
  theme(legend.position = "top") +
  facet_wrap(~Model)



# Assuming 'predictions' and 'actuals' are data frames containing the model predictions and actual values
predictions <- data.frame(Value = c(pred_tmb1_test, pred_tmb2_test, pred_tmb3_test), 
                          Type = "Predicted", 
                          Model = rep(c("tmb1", "tmb2", "tmb3"), each = nrow(test_data)))
actuals <- data.frame(Value = rep(test_data$log_return, 3), 
                      Type = "Actual", 
                      Model = rep(c("tmb1", "tmb2", "tmb3"), each = nrow(test_data)))

predicted_vs_actual <- rbind(predictions, actuals)

ggplot(predicted_vs_actual, aes(x = Value, fill = Type)) +
  geom_density(alpha = 0.7) +
  facet_wrap(~Model) +
  labs(title = "Predicted vs Actual Values", x = "Values", y = "Density") +
  theme_minimal() +
  theme(legend.position = "top")





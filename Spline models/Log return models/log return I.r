source("requirements.R")

plan(sequential)
gc()
options(future.seed = TRUE, future.rng.onMisuse="ignore")

getSymbols("FTSE", from = "2014-01-01")

FTSE_df <- as.data.frame(`FTSE`)
FTSE_df <- na.aggregate(FTSE_df)
colnames(FTSE_df) <- c("open", "high", "low", "close", "volume", "adjusted")
FTSE_df$time <- 1:nrow(FTSE_df)
FTSE_df$date <- index(`FTSE`)
FTSE_df$year <- format(FTSE_df$date, "%Y")
FTSE_df$month <- format(FTSE_df$date, "%m")
FTSE_df$log_return <- c(NA, diff(log(FTSE_df$close)))

fit_data <- na.omit(FTSE_df)
fit_data$close_open_ratio <- fit_data$close / fit_data$open
fit_data$volume <- lag(fit_data$volume, 1)
fit_data$close_open_ratio <- lag(fit_data$close_open_ratio, 1)
fit_data$close_lag <- lag(fit_data$close)
fit_data$log_ret_lag <- lag(fit_data$log_return)

calculate_rmse <- function(actual, predicted) {
  sqrt(mean((actual - predicted)^2))
}

set.seed(123)
fit_data <- fit_data[order(fit_data$time), ]
total_rows <- nrow(fit_data)
train_start_size <- floor(0.7 * total_rows) # Initial size of the training set
validation_set_size <- 100 # Fixed size for each validation set
n_folds <- floor((total_rows - train_start_size) / validation_set_size)
# Initialize RMSE storage
Linear_rmse_values <- c()
GAM_rmse_values <- c()

calculate_features_training <- function(data) {
  ts_data <- ts(data$log_return, frequency = 252) # Adjust frequency as needed
  decomposed <- stl(ts_data, s.window = "periodic", robust = TRUE)
  data$seasonal <- decomposed$time.series[, "seasonal"]
  data$trend <- decomposed$time.series[, "trend"]
  data$ema_20 <- SMA(data$close, n = 5)
  data$rsi_20 <- RSI(data$close, n = 5)
  data$hist_vol_20 <- runSD(data$log_return, n = 5, sample = FALSE)
  return(data)
}

initialize_validation_features <- function(train_data, valid_data) {
  last_train_values <- tail(train_data, 1)
  # Assume trend and seasonal components are constant over the validation period
  valid_data$seasonal <- rep(last_train_values$seasonal, nrow(valid_data))
  valid_data$trend <- rep(last_train_values$trend, nrow(valid_data))
  valid_data$ema_20 <- SMA(c(last_train_values$ema_20, valid_data$close), n = 5)[-1]
  valid_data$rsi_20 <- RSI(c(last_train_values$rsi_20, valid_data$close), n = 5)[-1]
  valid_data$hist_vol_20 <- runSD(c(last_train_values$hist_vol_20, valid_data$log_return), n = 5, sample = FALSE)[-1]
  return(valid_data)
}

# Initialize RMSE storage
Linear_rmse_values <- c()
GAM_rmse_values <- c()
fold_data_list <- list()

# Determine the number of cores to use
n_cores <- parallel::detectCores(logical = TRUE) - 7
plan(multisession, workers = n_cores)

handlers(global = TRUE)
handler_progress()

pb <- progress_bar$new(
  format = "  Folding [:bar] :percent :etas",
  total = n_folds, clear = FALSE, width = 30
)

process_fold <- function(i) {
  current_train_end <- min(train_start_size + (i - 1) * validation_set_size, total_rows - validation_set_size)
  current_train_data <- fit_data[1:current_train_end, ]
  
  current_train_data <- calculate_features_training(current_train_data)
  current_valid_data <- fit_data[(current_train_end + 1):(current_train_end + validation_set_size), ]
  current_valid_data <- initialize_validation_features(current_train_data, current_valid_data)
  
  current_train_data <- na.omit(current_train_data)
  current_valid_data <- na.omit(current_valid_data)
  
  pb$tick()
  
  Linear_cv <- lm(log_return ~ time + volume + close_open_ratio + 
                  ema_20 + rsi_20 + hist_vol_20 +
                  trend + seasonal, 
                  data = current_train_data)
  
  GAM_cv <- glmmTMB(log_return ~ s(time) + s(volume) + s(close_open_ratio) + 
                    s(ema_20) + s(hist_vol_20) + s(rsi_20) +
                    s(trend) + seasonal,
                    disp =~ month,
                    data = current_train_data, 
                    family = gaussian(link = "identity"),
                    REML = TRUE)
  
  Linear_predictions <- predict(Linear_cv, current_valid_data, type="response")
  GAM_predictions <- predict(GAM_cv, current_valid_data, type="response", allow.new.levels = TRUE)
  
  Linear_rmse_values <- c(Linear_rmse_values, calculate_rmse(current_valid_data$log_return, Linear_predictions))
  GAM_rmse_values <- c(GAM_rmse_values, calculate_rmse(current_valid_data$log_return, GAM_predictions))
  
  cat(sprintf("Fold %d: Linear Model RMSE = %f, GAM Model RMSE = %f\n", i, tail(Linear_rmse_values, 1), tail(GAM_rmse_values, 1)))
  
  Linear_rmse_value <- calculate_rmse(current_valid_data$log_return, Linear_predictions)
  GAM_rmse_value <- calculate_rmse(current_valid_data$log_return, GAM_predictions)
  
  # Inefficient but works for parallelized CV loop 
  return(list(
    fold = i,
    valid_data = current_valid_data,
    Linear_predictions = Linear_predictions,
    GAM_predictions = GAM_predictions,
    Linear_rmse = Linear_rmse_value,
    GAM_rmse = GAM_rmse_value
  ))
}
  
results <- with_progress({
  future_map_dfr(1:n_folds, process_fold)
})

average_Linear_rmse <- mean(results$Linear_rmse)
average_GAM_rmse <- mean(results$GAM_rmse)

cat(sprintf("Average Linear Model RMSE: %f\n", average_Linear_rmse))
cat(sprintf("Average GAM Model RMSE: %f\n", average_GAM_rmse))

results_long <- tidyr::pivot_longer(results, cols = c(Linear_rmse, GAM_rmse), names_to = "Model", values_to = "RMSE")

ggplot(results_long, aes(x = as.factor(fold), y = RMSE, color = Model)) +
  geom_point() +
  geom_line(aes(group = Model), position = position_dodge(width = 0.25)) +
  labs(title = "RMSE Across Folds for Linear and GAM Models", x = "Fold", y = "RMSE") +
  theme_minimal() +
  scale_color_manual(values = c("Linear_rmse" = "darkblue", "GAM_rmse" = "darkorange"))

fold_predictions_actuals <- data.frame(
   fold = results_long$fold, # Fold numbers
   actual_log_return = results_long$valid_data$log_return, # Actual log returns
   linear_prediction = results_long$Linear_predictions, # Predictions from the Linear model
   gam_prediction = results_long$GAM_predictions # Predictions from the GAM model
)


fold_predictions_actuals_long <- fold_predictions_actuals %>%
  tidyr::pivot_longer(
    cols = c(linear_prediction, gam_prediction),
    names_to = "Model",
    values_to = "Predicted"
  ) %>%
  mutate(Model = ifelse(Model == "linear_prediction", "Linear", "GAM"))

ggplot(fold_predictions_actuals_long, aes(x = actual_log_return, y = Predicted, color = Model)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE) +
  labs(title = "Predicted vs. Actual Log Returns with LOESS", x = "Actual Log Return", y = "Predicted Log Return") +
  theme_minimal() +
  scale_color_manual(values = c("Linear" = "darkblue", "GAM" = "darkorange"))

ggplot(fold_predictions_actuals_long, aes(x = Predicted, fill = Model)) +
  geom_density(alpha = 0.6) +
  labs(title = "Density of Predicted Log Returns for Each Model", x = "Predicted Log Return", y = "Density") +
  theme_minimal() +
  scale_fill_manual(values = c("Linear" = "darkblue", "GAM" = "darkorange"))


ggplot(fold_predictions_actuals_long, aes(x = actual_log_return, y = Predicted, color = Model)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Predicted vs. Actual Log Returns with LOESS", x = "Actual Log Return", y = "Predicted Log Return") +
  theme_minimal() +
  scale_color_manual(values = c("Linear" = "darkblue", "GAM" = "darkorange"))


ggplot(fold_predictions_actuals_long, aes(x = Predicted, fill = Model)) +
  geom_density(alpha = 0.6) +
  # Add actual densities
  geom_density(data = fold_predictions_actuals_long, aes(x = actual_log_return, y = ..density..), color = "black", fill = NA) +
  labs(title = "Density of Predicted Log Returns vs. Actual", x = "Log Return", y = "Density") +
  theme_minimal() +
  scale_fill_manual(values = c("Linear" = "darkblue", "GAM" = "darkorange"))


source("requirements.R")

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

covariates <- c("time", "N100_RSI_lag", "N100_log_return_lag", "N100_trend", "N225_RSI", "N225_log_return")

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

calculate_n_knots <- function(smooth_term) {
  if (!is.null(smooth_term[["rank"]])) {
    # For tensor product smooths, use the first element of 'rank'
    return(smooth_term[["rank"]][1])
  } else {
    return(smooth_term[["bs.dim"]] - smooth_term[["null.space.dim"]])
  }
}

generate_smooth_terms <- function(covariate, data, bs = "cs", k) {
  # Ensure that k is provided
  if (is.null(k)) stop("The 'k' argument (number of knots) must be provided.")
  formula_str <- paste0("s(", covariate, ", bs='", bs, "', k=", k, ")")
  sm_formula <- eval(parse(text = formula_str))
  sm <- mgcv::smoothCon(sm_formula, data = data, absorb.cons=TRUE)[[1]]
  re <- mgcv::smooth2random(sm, "", type = 2)
  pred_matrix <- s2rPred(sm, re, data = data)
  
  # Return the necessary components
  list(Xf = pred_matrix$Xf, Xr = pred_matrix$rand[[1]], n_knots = calculate_n_knots(sm))
}
results <- future_map(covariates, ~generate_smooth_terms(.x, data = merged_df, bs = "cs", k = 5))

for (i in seq_along(covariates)) {
  covariate <- covariates[i]
  merged_df[[paste0("Xf_", covariate)]] <- results[[i]]$Xf
  merged_df[[paste0("Xr_", covariate)]] <- results[[i]]$Xr
  merged_df[[paste0("n_knots_", covariate)]] <- results[[i]]$n_knots
}
n_knots_list <- setNames(lapply(results, `[[`, "n_knots"), covariates)

augment_design_matrix <- function(Xr, lambda) {
  if (is.null(Xr) || nrow(Xr) == 0 || ncol(Xr) == 0) {
    stop("Input matrix Xr is NULL or empty")
  }
  tryCatch({
    augmented_matrix <- rbind(Xr, sqrt(lambda) * diag(ncol(Xr)))
    if (is.null(augmented_matrix) || nrow(augmented_matrix) == 0) {
      stop("Augmented matrix is NULL or empty after applying ridge penalty")
    }
  }, error = function(e) {
    stop(paste("Error in augment_design_matrix: ", e$message))
  })
  
  return(augmented_matrix)
}

augment_response <- function(y, augmented_length) {
  c(y, rep(0, augmented_length - length(y)))
}


augment_data <- function(data, lambda, n_knots_list) {
  augmented_matrices <- list()
  
  for (name in names(n_knots_list)) {
    xr_var <- data[[paste0("Xr_", name)]]
    if (is.null(xr_var)) {
      stop(paste0("Variable ", paste0("Xr_", name), " is NULL"))
    }
    augmented_matrix <- augment_design_matrix(xr_var, lambda)
    augmented_matrices[[name]] <- augmented_matrix
  }
  
  y_augmented <- augment_response(data$N100_log_return, max(sapply(augmented_matrices, nrow)))
  
  augmented_data <- data.frame(N100_log_return = y_augmented)  # Ensure correct response variable name
  for (name in names(n_knots_list)) {
    augmented_matrix <- augmented_matrices[[name]]
    augmented_data[[paste0("Xr_", name)]] <- I(augmented_matrix)
  }
  return(augmented_data)
}

gcv_lambda <- function(lambda, n_knots_list, data) {
  augmented_data <- augment_data(data, lambda, n_knots_list)
  
  gcv_model <- glmmTMB(formula = N100_log_return ~ 1 + Xr_time + 
                              Xr_N100_RSI_lag + Xr_N225_log_return, 
                       data = augmented_data, family = gaussian(link="identity"))
  
  log_likelihood <- logLik(gcv_model)
  edf <- attr(log_likelihood, "df")
  deviance_val <- deviance(gcv_model)
  N <- nrow(augmented_data)
  gcv_score <- deviance_val / ((1 - edf/N)^2)
  
  return(gcv_score)
}

cross_validate_lambda_with_early_stopping <- function(data, lambda_values, n_knots_list, early_stopping_rounds = 3) {
  require(future.apply)
  plan(multisession, workers = 5) 
  best_score <- Inf
  scores <- numeric(length(lambda_values))
  no_improvement_count <- 0
  
  lambda_scores <- future_map(lambda_values, function(lambda) {
    gcv_lambda(lambda, n_knots_list, data)
  })
  
  for (i in seq_along(lambda_scores)) {
    score <- lambda_scores[[i]]
    scores[i] <- score
    
    if (score < best_score) {
      best_score <- score
      no_improvement_count <- 0
    } else {
      no_improvement_count <- no_improvement_count + 1
      if (no_improvement_count >= early_stopping_rounds) {
        message(sprintf("Early stopping after %d iterations", i))
        break
      }
    }
  }
  
  scores <- scores[1:i]  # Truncate the scores vector to actual number of iterations
  lambda_values_truncated <- lambda_values[1:i]  # Truncate lambda values to match the scores length
  
  list(optimal_lambda = lambda_values[which.min(scores)], gcv_scores = scores, lambda_values = lambda_values_truncated)
}

plan(sequential)
gc()
options(future.seed = TRUE, future.rng.onMisuse="ignore")

lambda_values <- seq(0.00001, 1, by = 0.05)

# Time the cross-validation process
timing <- system.time({
  cv_results <- cross_validate_lambda_with_early_stopping(merged_df, lambda_values, n_knots_list, early_stopping_rounds = 5)
})

# Extract the results
optimal_lambda <- cv_results$optimal_lambda
gcv_values <- cv_results$gcv_scores
lambda_values_used <- cv_results$lambda_values

# Print the timing results
print(timing)

merged_df_augmented <- augment_data(merged_df, optimal_lambda, n_knots_list)

plan(sequential)
gc()
options(future.seed = TRUE, future.rng.onMisuse="ignore")
n_cores <- parallel::detectCores(logical = TRUE) - 7
plan(multisession, workers = n_cores)

manual_model_formula <- formula(N100_log_return ~ 1 + Xr_time + 
                                  Xr_N100_RSI_lag + Xr_N225_log_return)
manual_model_family <- gaussian(link = "identity")

auto_model_formula <- formula(N100_log_return ~ s(time) + s(N100_RSI_lag) + s(N225_log_return))
auto_model_family <- gaussian(link = "identity")

set.seed(10)
train_proportion <- 0.6
valid_proportion <- 0.2
total_rows <- nrow(merged_df_clean)
train_size <- floor(train_proportion * total_rows)
valid_size <- floor(valid_proportion * total_rows)

train_data <- merged_df_clean[1:train_size, ]
valid_data <- merged_df_clean[(train_size + 1):(train_size + valid_size), ]
test_data <- merged_df_clean[(train_size + valid_size + 1):total_rows, ]

train_data_augmented <- merged_df_augmented[1:train_size, ]
valid_data_augmented <- merged_df_augmented[(train_size + 1):(train_size + valid_size), ]
test_data_augmented <- merged_df_augmented[(train_size + valid_size + 1):total_rows, ]

calculate_metrics <- function(actual, predicted) {
  sqrt(mean((predicted - actual)^2))
}

fit_and_evaluate_model <- function(model_formula, train_data, valid_data, family, use_REML = FALSE) {
  timing <- system.time({
    model <- glmmTMB(model_formula, data = train_data, family = family, REML = use_REML)
    predictions <- predict(model, newdata = valid_data, re.form = NA)
  })
  metrics <- list(
    RMSE = calculate_metrics(valid_data$N100_log_return, predictions),
    TimeTaken = timing["elapsed"]
  )
  return(metrics)
}

library(furrr)
library(glmmTMB)
library(ggplot2)
library(reshape2)

calculate_metrics <- function(actual, predicted) {
  sqrt(mean((predicted - actual)^2))
}

fit_and_evaluate_model <- function(model_formula, train_data, valid_data, family, use_REML = FALSE) {
  timing <- system.time({
    model <- glmmTMB(model_formula, data = train_data, family = family, REML = use_REML)
    predictions <- predict(model, newdata = valid_data, re.form = NA)
  })
  metrics <- list(
    RMSE = calculate_metrics(valid_data$N100_log_return, predictions),
    TimeTaken = timing["elapsed"]
  )
  return(metrics)
}

cross_validate_models <- function(data, augmented_data, folds) {
  metrics_list <- list()
  train_sets <- list()
  valid_sets <- list()
  
  for (fold in unique(folds)) {
    validation_indexes <- which(folds == fold)
    cv_train <- data[-validation_indexes, ]
    cv_valid <- data[validation_indexes, ]
    cv_train_aug <- augmented_data[-validation_indexes, ]
    cv_valid_aug <- augmented_data[validation_indexes, ]
    
    manual_model_metrics <- fit_and_evaluate_model(manual_model_formula, cv_train_aug, cv_valid_aug, manual_model_family)
    auto_model_metrics <- fit_and_evaluate_model(auto_model_formula, cv_train, cv_valid, auto_model_family, use_REML = FALSE)
    
    # Print RMSE and time taken for each fold
    cat(sprintf("Fold %d - ManualModel RMSE: %f, Time Taken (s): %f\n", fold, manual_model_metrics$RMSE, manual_model_metrics$TimeTaken))
    cat(sprintf("Fold %d - AutoModel RMSE: %f, Time Taken (s): %f\n", fold, auto_model_metrics$RMSE, auto_model_metrics$TimeTaken))
    
    metrics_list[[fold]] <- list(
      ManualModel = manual_model_metrics,
      AutoModel = auto_model_metrics
    )
    
    # Store the train and validation sets
    train_sets[[fold]] <- list(
      cv_train = cv_train,
      cv_train_aug = cv_train_aug
    )
    valid_sets[[fold]] <- list(
      cv_valid = cv_valid,
      cv_valid_aug = cv_valid_aug
    )
  }
  
  return(list(metrics = metrics_list, train_sets = train_sets, valid_sets = valid_sets))
}

k <- 5
folds <- cut(seq(1, nrow(merged_df)), breaks = k, labels = FALSE)

# Time the cross-validation process
total_timing <- system.time({
  cv_results <- cross_validate_models(merged_df, merged_df_augmented, folds)
})

# Print the total timing results
cat(sprintf("Total cross-validation time: %f seconds\n", total_timing["elapsed"]))

# Combine the train and validation sets into "full_train_data" sets
full_train_data <- do.call(rbind, lapply(cv_results$train_sets, function(x) x$cv_train))
full_train_data_augmented <- do.call(rbind, lapply(cv_results$train_sets, function(x) x$cv_train_aug))

# Define the test data (20% of the original data)
set.seed(123)  # For reproducibility
test_indexes <- sample(1:nrow(merged_df), size = floor(0.2 * nrow(merged_df)))
test_data <- merged_df[test_indexes, ]
test_data_augmented <- merged_df_augmented[test_indexes, ]

# Ensure the training data does not include the test data
full_train_data <- full_train_data[!rownames(full_train_data) %in% test_indexes, ]
full_train_data_augmented <- full_train_data_augmented[!rownames(full_train_data_augmented) %in% test_indexes, ]

# Print the dimensions to verify
cat(sprintf("Full Train Data Dimensions: %d x %d\n", nrow(full_train_data), ncol(full_train_data)))
cat(sprintf("Full Train Data Augmented Dimensions: %d x %d\n", nrow(full_train_data_augmented), ncol(full_train_data_augmented)))
cat(sprintf("Test Data Dimensions: %d x %d\n", nrow(test_data), ncol(test_data)))
cat(sprintf("Test Data Augmented Dimensions: %d x %d\n", nrow(test_data_augmented), ncol(test_data_augmented)))

fit_and_predict <- function(model_formula, family, train_data, valid_data, augmented = FALSE, augmented_train_data = NULL, use_REML = FALSE) {
  data_to_use <- if (augmented) augmented_train_data else train_data
  
  # Time the model fitting
  training_time <- system.time({
    model <- glmmTMB(model_formula, data = data_to_use, family = family, REML = use_REML)
  })
  
  # Time the predictions
  prediction_time <- system.time({
    predictions <- predict(model, newdata = valid_data, re.form = NA)
  })
  
  list(
    Predictions = predictions,
    TrainingTime = training_time["elapsed"],
    PredictionTime = prediction_time["elapsed"]
  )
}

# Train and evaluate the final models
final_manual_model <- fit_and_predict(manual_model_formula, manual_model_family, full_train_data_augmented, test_data_augmented, augmented = TRUE, augmented_train_data = full_train_data_augmented, use_REML = FALSE)
final_auto_model <- fit_and_predict(auto_model_formula, auto_model_family, full_train_data, test_data, augmented = FALSE, use_REML = FALSE)

# Print final model performance
cat(sprintf("Final Manual Model - RMSE: %f, Training Time (s): %f, Prediction Time (s): %f\n", 
            calculate_metrics(test_data$N100_log_return, final_manual_model$Predictions), 
            final_manual_model$TrainingTime, 
            final_manual_model$PredictionTime))
cat(sprintf("Final Auto Model - RMSE: %f, Training Time (s): %f, Prediction Time (s): %f\n", 
            calculate_metrics(test_data$N100_log_return, final_auto_model$Predictions), 
            final_auto_model$TrainingTime, 
            final_auto_model$PredictionTime))


# Combine predictions into a list
model_predictions <- list(
  ManualModel = final_manual_model,
  AutoModel = final_auto_model
)

plan(sequential)
gc()
options(future.seed = TRUE, future.rng.onMisuse="ignore")
plan(multisession, workers = 2)

plot_predictions <- function(model_predictions, test_data, models_to_train) {
  prediction_df <- data.frame(Actual = test_data$N100_log_return)
  rmse_values <- list()
  timing_values <- list()
  
  for (model_name in names(models_to_train)) {
    predictions <- model_predictions[[model_name]]$Predictions
    prediction_df[[model_name]] <- predictions
    
    # Calculate RMSE for the current model's predictions
    rmse_value <- calculate_metrics(prediction_df$Actual, predictions)
    rmse_values[[model_name]] <- rmse_value
    
    # Capture timing information
    training_time <- model_predictions[[model_name]]$TrainingTime
    prediction_time <- model_predictions[[model_name]]$PredictionTime
    timing_values[[model_name]] <- list(TrainingTime = training_time, PredictionTime = prediction_time)
    
    # Print RMSE value and timing information in the console
    cat(sprintf("RMSE for %s: %f\n", model_name, rmse_value))
    cat(sprintf("Training time for %s: %f seconds\n", model_name, training_time))
    cat(sprintf("Prediction time for %s: %f seconds\n", model_name, prediction_time))
  }
  
  # Update model names for the plot
  names(rmse_values)[names(rmse_values) == "ManualModel"] <- "Ridge"
  names(rmse_values)[names(rmse_values) == "AutoModel"] <- "ISSD"
  
  # Reshape the data for plotting
  prediction_df_long <- melt(prediction_df, id.vars = "Actual", variable.name = "Model", value.name = "Predicted")
  prediction_df_long$Model <- factor(prediction_df_long$Model, levels = names(models_to_train))
  levels(prediction_df_long$Model) <- c("Ridge", "ISSD") # Rename levels to match desired output
  
  # Create the plot
  p <- ggplot(prediction_df_long, aes(x = Actual, y = Predicted, color = Model, shape = Model)) +
    geom_point(position = position_jitter(width = 0.05, height = 0), alpha = 0.4, size = 1) +
    geom_smooth(method = "loess", se = TRUE, linewidth = 1, aes(color = Model)) +
    geom_abline(linetype = "dashed", color = "red", linewidth = 1) +
    scale_shape_manual(values = c(16, 17)) +
    scale_color_manual(values = c("Ridge" = "darkblue", "ISSD" = "darkorange"), labels = c("Ridge", "ISSD")) +
    theme_minimal(base_size = 12) +
    labs(
      x = "Actual", 
      y = "Predicted", 
      color = "Model",
      title = "Model Predictions vs Actual Values on Test Data",
      subtitle = "Smoothed Regression Line with Confidence Interval and Theoretical Line"
    ) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 14),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.position = "bottom"
    ) +
    guides(
      color = guide_legend(override.aes = list(alpha = 1, size = 3)),
      shape = guide_legend(override.aes = list(size = 3))
    )
  
  # Add RMSE values as annotations on the plot
  num_models <- length(rmse_values)
  model_names <- names(rmse_values)
  for (i in seq_along(model_names)) {
    model_name <- model_names[i]
    
    vertical_offset <- 3.7 - (i-1) * 1.7
    
    p <- p + annotate("text", x = Inf, y = Inf, 
                      label = paste(model_name, "RMSE:", round(rmse_values[[model_name]], 3)), 
                      hjust = 3, # Adjust if necessary to move text left from the right edge
                      vjust = vertical_offset, # Dynamic vertical adjustment based on model index
                      color = ifelse(model_name == "Ridge", "darkblue", "darkorange"), size = 10)
  }
  
  print(p)
}

# Plot the predictions
plot_predictions(model_predictions, test_data, list(ManualModel = manual_model_formula, AutoModel = auto_model_formula))


plot_predictions(model_predictions, test_data, models_to_train)



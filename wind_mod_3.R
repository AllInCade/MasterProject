library(parallel)
library(doParallel)
library(dplyr)
library(zoo)
library(lubridate)
library(readr)
library(mgcv)
library(glmmTMB)
library(randomForest)
#library(statmod)
#library(tweedie)
library(forecast)
library(ggplot2)
library(DHARMa)
library(reshape2)
library(purrr)
library(furrr)
library(future)
library(future.apply)

plan(sequential)
gc()
options(future.seed = TRUE, future.rng.onMisuse="ignore")
gc()
wind_data <- read_csv("SB_Aeolus_20230519.csv")
wind_data <- wind_data %>%
  mutate(datetime_UTC = as.POSIXct(datetime_UTC, format="%Y-%m-%d %H:%M:%S")) %>%
  filter(datetime_UTC >= as.POSIXct("2007-09-01") & datetime_UTC <= as.POSIXct("2007-10-18")) %>%
  filter(!(tolower(abb) == c("sm", "46029") & tolower(name) == "santa_maria")) %>%
  mutate(numeric_time = as.numeric(difftime(datetime_UTC, as.POSIXct("1970-01-01"), units = "hours")),
         time = numeric_time - min(numeric_time)) %>%
  mutate(across(everything(), forecast::na.interp)) %>%
  filter(spd >= 0.1) %>%
  mutate(spd_rolling_avg = rollapply(spd, width = 7, FUN = function(x) mean(x, na.rm = TRUE), 
                                     align = 'right', fill = NA, partial = TRUE))

wind_data <- wind_data[order(wind_data$datetime_UTC), ]
date_row <- min(which(wind_data$datetime_UTC == as.POSIXct("2007-10-12")))

rm(wind_data)
gc()
sample_interval <- 1  # Can be adjusted as needed
wind_data_clean <- wind_data_clean[seq(1, nrow(wind_data_clean), sample_interval), ]

wind_data_clean$datetime_UTC <- as.POSIXct(wind_data_clean$datetime_UTC)
wind_ts <- ts(wind_data_clean$spd, frequency = 48) # capture diurnal pattern

no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)

fourier_terms <- 1:5
models <- list()
models <- foreach(k = fourier_terms, .packages = 'forecast') %dopar% {
  auto.arima(wind_ts, seasonal = TRUE, stepwise = FALSE, approximation = FALSE,
             D = 1,  # Level of differencing
             max.P = 2, max.Q = 2,  # Adjust based on the length of the season
             xreg = fourier(wind_ts, K=k))  # Fourier terms for complex seasonality
}
stopCluster(cl)

best_model <- models[[which.min(sapply(models, AIC))]]
checkresiduals(best_model)
residuals_data <- residuals(best_model)
wind_data_clean$residuals_sarima <- residuals_data

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

covariates <- c("residuals_sarima", "spd_rolling_avg", "Vtau", "Utau", "U", "V")
generate_smooth_terms <- function(covariate, data) {
  formula_str <- paste0("s(", covariate, ")")
  sm_formula <- eval(parse(text = formula_str))
  sm <- mgcv::smoothCon(sm_formula, absorb.cons=TRUE, data = data)[[1]]
  re <- mgcv::smooth2random(sm, "", type = 2)
  pred_matrix <- s2rPred(sm, re, data = data)
  list(Xf = pred_matrix$Xf, Xr = pred_matrix$rand[[1]], n_knots = calculate_n_knots(sm))
}
results <- future_map(covariates, ~generate_smooth_terms(.x, data = wind_data_clean)) 

for (i in seq_along(covariates)) {
  covariate <- covariates[i]
  wind_data_clean[[paste0("Xf_", covariate)]] <- results[[i]]$Xf
  wind_data_clean[[paste0("Xr_", covariate)]] <- results[[i]]$Xr
  wind_data_clean[[paste0("n_knots_", covariate)]] <- results[[i]]$n_knots
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
  
  y_augmented <- augment_response(data$spd, max(sapply(augmented_matrices, nrow)))
  
  augmented_data <- data.frame(spd = y_augmented)  # Ensure correct response variable name
  for (name in names(n_knots_list)) {
    augmented_matrix <- augmented_matrices[[name]]
    augmented_data[[paste0("Xr_", name)]] <- I(augmented_matrix)
  }
  return(augmented_data)
}

gcv_lambda <- function(lambda, n_knots_list, data) {
  augmented_data <- augment_data(data, lambda, n_knots_list)
  
  gcv_model <- glmmTMB(formula = spd ~ 1 + Xr_Utau + Xr_Vtau + Xr_spd_rolling_avg +
                         + Xr_residuals_sarima + Xr_U + Xr_V, 
                       data = augmented_data, family = gaussian(link="identity"))
  
  log_likelihood <- logLik(gcv_model)
  edf <- attr(log_likelihood, "df")
  deviance_val <- deviance(gcv_model)
  N <- nrow(augmented_data)
  gcv_score <- deviance_val / ((1 - edf/N)^2)
  
  return(gcv_score)
}

cross_validate_lambda <- function(data, lambda_values, n_knots_list) {
  plan(multisession, workers = 5)  # Setting up parallel processing
  gcv_scores <- future_sapply(lambda_values, function(lambda) {
    gcv_lambda(lambda, n_knots_list, data)
  })
  
  list(optimal_lambda = lambda_values[which.min(gcv_scores)], gcv_scores = gcv_scores)
}

plan(sequential)
gc()
options(future.seed = TRUE, future.rng.onMisuse="ignore")
plan(multisession, workers = 5)

lambda_values <- seq(0.001, 0.2, by = 0.005)
cv_results <- cross_validate_lambda(wind_data_clean, lambda_values, n_knots_list)
optimal_lambda <- cv_results$optimal_lambda
gcv_values <- cv_results$gcv_scores
wind_data_clean_augmented <- augment_data(wind_data_clean, optimal_lambda, n_knots_list)

lambda_performance_df <- data.frame(lambda = lambda_values, GCV_Score = gcv_values)
ggplot(lambda_performance_df, aes(x = lambda, y = GCV_Score)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Lambda vs. GCV Score", x = "Lambda", y = "GCV Score")


plan(sequential)
gc()
options(future.seed = TRUE, future.rng.onMisuse="ignore")
plan(multisession, workers = 2)  # This will choose multicore on Unix/Linux and multisession on Windows

manual_model_formula <- formula(spd ~ 1 + Xr_Utau + Xr_Vtau + Xr_spd_rolling_avg +
                                  + Xr_residuals_sarima + Xr_U + Xr_V)
manual_model_family <- gaussian(link = "identity")

auto_model_formula <- formula(spd ~ s(Utau) + s(Vtau) + s(spd_rolling_avg) + s(residuals_sarima) + s(U) + s(V))
auto_model_family <- gaussian(link = "identity")

set.seed(10)
train_proportion <- 0.6
valid_proportion <- 0.2
total_rows <- nrow(wind_data)
train_size <- floor(train_proportion * total_rows)
valid_size <- floor(valid_proportion * total_rows)

if (train_size + valid_size >= date_row) {
  train_size <- date_row - 1 - floor(valid_proportion * total_rows) # might need tuning
}

train_data <- wind_data[1:train_size, ]
valid_data <- wind_data[(train_size + 1):(train_size + valid_size), ]
test_data <- wind_data[(train_size + valid_size + 1):total_rows, ]

train_data_augmented <- wind_data_clean_augmented[1:train_size, ]
valid_data_augmented <- wind_data_clean_augmented[(train_size + 1):(train_size + valid_size), ]
test_data_augmented <- wind_data_clean_augmented[(train_size + valid_size + 1):total_rows, ]

calculate_metrics <- function(actual, predicted) {
  sqrt(mean((predicted - actual)^2))
}

fit_and_evaluate_model <- function(model_formula, train_data, valid_data, family) {
  start_time <- Sys.time()
  model <- glmmTMB(model_formula, data = train_data, family = family)
  predictions <- predict(model, newdata = valid_data, re.form = NA)
  end_time <- Sys.time()
  metrics <- list(RMSE = calculate_metrics(valid_data$spd, predictions),
                  TimeTaken = end_time - start_time)
  return(metrics)
}

cross_validate_models <- function(data, augmented_data, folds) {
  metrics_list <- future_map(unique(folds), function(fold) {
    validation_indexes <- which(folds == fold)
    cv_train <- data[-validation_indexes, ]
    cv_valid <- data[validation_indexes, ]
    cv_train_aug <- augmented_data[-validation_indexes, ]
    cv_valid_aug <- augmented_data[validation_indexes, ]
    list(
      ManualModel = fit_and_evaluate_model(manual_model_formula, cv_train_aug, cv_valid_aug, manual_model_family),
      AutoModel = fit_and_evaluate_model(auto_model_formula, cv_train, cv_valid, auto_model_family)
    )
  })
  return(metrics_list)
}

k <- 5
folds <- cut(seq(1, nrow(wind_data_clean)), breaks = k, labels = FALSE)
results <- cross_validate_models(wind_data_clean, wind_data_clean_augmented, folds)
avg_metrics <- map(results, ~sapply(.x, function(x) sapply(x, mean)))
print(avg_metrics)

proportions <- seq(0.1, 1, by = 0.1)
fit_and_evaluate_on_subset <- function(model_formula, model_family, full_train_data, valid_data, proportion, augmented = FALSE, augmented_train_data = NULL) {
  set.seed(123)  # Ensure reproducibility
  sample_indices <- sample(nrow(full_train_data), size = floor(proportion * nrow(full_train_data)))
  train_subset <- if (augmented && !is.null(augmented_train_data)) {
    augmented_train_data[sample_indices, ]
  } else {
    full_train_data[sample_indices, ]
  }
  fit_and_evaluate_model(model_formula, train_subset, valid_data, model_family)
}

learning_curve_data <- function(data, valid_data, proportions, augmented_data = NULL) {
  future_map(proportions, function(prop) {
    list(
      ManualModel = fit_and_evaluate_on_subset(manual_model_formula, manual_model_family, data, valid_data, prop, augmented = TRUE, augmented_data),
      AutoModel = fit_and_evaluate_on_subset(auto_model_formula, auto_model_family, data, valid_data, prop)
    )
  })
}

plan(sequential)
gc()
options(future.seed = TRUE, future.rng.onMisuse="ignore")
plan(multisession, workers = 4)

plot_learning_curve <- function(learning_curve_data, proportions) {
  learning_curve_df <- data.frame(
    Proportion = rep(proportions, each = 2),
    Model = rep(c("ManualModel", "AutoModel"), times = length(proportions)),
    RMSE = unlist(lapply(learning_curve_data, function(x) c(x$ManualModel$RMSE, x$AutoModel$RMSE)))
  )
  
  ggplot(learning_curve_df, aes(x = Proportion, y = RMSE, color = Model)) +
    geom_point() +
    geom_smooth(method = "loess", se = FALSE) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    ) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    labs(
      title = "Learning Curve",
      subtitle = "Performance of Manual and Automatic Models with Increasing Training Data Proportion",
      x = "Proportion of Training Data",
      y = "RMSE",
      color = "Model"
    )
}

lc_data <- learning_curve_data(wind_data_clean, valid_data, proportions, wind_data_clean_augmented)
plot_learning_curve(lc_data, proportions)

fit_and_predict <- function(model_formula, model_family, train_data, test_data, augmented = FALSE, augmented_train_data = NULL) {
  model_data <- if (augmented && !is.null(augmented_train_data)) augmented_train_data else train_data
  model <- glmmTMB(model_formula, data = model_data, family = model_family)
  predictions <- predict(model, newdata = test_data, re.form = NA)
  return(predictions)
}

full_train_data <- rbind(train_data, valid_data)
full_train_data_augmented <- rbind(train_data_augmented, valid_data_augmented)

plan(sequential)
gc()
options(future.seed = TRUE, future.rng.onMisuse="ignore")# If using augmented data
plan(multisession, workers = 2)
models_to_train <- list(
  ManualModel = list(formula = manual_model_formula, family = manual_model_family, augmented = TRUE),
  AutoModel = list(formula = auto_model_formula, family = auto_model_family, augmented = FALSE)
)

model_predictions <- future_map(models_to_train, function(model_info) {
  fit_and_predict(model_info$formula, model_info$family, full_train_data, test_data, augmented = model_info$augmented, full_train_data_augmented)
})

plan(sequential)
gc()
options(future.seed = TRUE, future.rng.onMisuse="ignore")# If using augmented data
plan(multisession, workers = 2)

plot_predictions <- function(model_predictions, test_data, models_to_train) {
  prediction_df <- data.frame(Actual = test_data$spd)
  
  for (i in seq_along(models_to_train)) {
    model_name <- names(models_to_train)[i]
    prediction_df[[model_name]] <- model_predictions[[i]]
  }
  
  prediction_df_long <- reshape2::melt(prediction_df, id.vars = "Actual", variable.name = "Model", value.name = "Predicted")
  ggplot(prediction_df_long, aes(x = Actual, y = Predicted, color = Model, shape = Model)) +
    geom_point(position = position_jitter(width = 0.05, height = 0), alpha = 0.2, size = 1) +
    geom_smooth(method = "loess", se = TRUE, linewidth = 1, color = "black") +  # Enhanced LOESS curve
    geom_abline(linetype = "dashed", color = "red", linewidth = 1) +  # Theoretical line x = y
    scale_shape_manual(values = c(16, 17)) +  # Different shapes for different models
    scale_color_manual(values = c("darkblue", "darkorange")) +  # Distinct color palette
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
      color = guide_legend(override.aes = list(alpha = 1, size = 2)),
      shape = guide_legend(override.aes = list(size = 2))
    )
}

plot_predictions(model_predictions, test_data, models_to_train)


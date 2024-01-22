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

temp_data <- read_csv("weather-anomalies-1964-2013.csv")
temp_data$date <- as.Date(temp_data$date_str, format="%Y-%m-%d")
temp_data$year <- as.integer(format(temp_data$date, "%Y"))

us_latitude_bounds <- c(24.396308, 49.384358)
us_longitude_bounds <- c(-125.001650, -66.934570)

temp_data_filtered <- temp_data %>%
  filter(year >= 2000,
         latitude >= us_latitude_bounds[1] & latitude <= us_latitude_bounds[2],
         longitude >= us_longitude_bounds[1] & longitude <= us_longitude_bounds[2])

sample_size = 1000

temp_data_stratified <- temp_data_filtered %>%
  group_by(year, type) %>%
  slice_sample(n = sample_size, replace = TRUE) %>%
  arrange(date)  


extreme_anomalies <- temp_data_stratified %>%
  filter(abs(degrees_from_mean) > quantile(abs(degrees_from_mean), 0.90))

ggplot(extreme_anomalies, aes(x = longitude, y = latitude)) +
  geom_point(aes(color = degrees_from_mean), shape = 3, size = 2, alpha = 0.8) +
  scale_color_gradientn(colors = c("blue", "green", "yellow", "red"), 
                        limits = c(-55, 55), na.value = "grey50") +
  labs(title = "Geographic Distribution of 90th Percentile Temperature Anomalies in the Contiguous US",
       subtitle = "Filtered Data from 2000 Onwards",
       x = "Longitude", y = "Latitude", color = "Temp Anomaly (°C)") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12),
    legend.position = "bottom"
  ) +
  coord_fixed(ratio = 1)

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

# Parallelization setup
plan(sequential)
gc()
options(future.seed = TRUE, future.rng.onMisuse="ignore")
plan(multisession, workers = 5)


rm(temp_data)
temp_data <- temp_data_stratified
colnames(temp_data)[11] <- "time"
temp_data$time <- as.numeric(temp_data$time - 10956)
#write_csv(temp_data, "tempdata.csv")

covariates <- list("time", "min_temp", "max_temp", c("latitude", "longitude"))

calculate_n_knots <- function(smooth_term) {
  if (!is.null(smooth_term[["rank"]])) {
    # For tensor product smooths, use the first element of 'rank'
    return(smooth_term[["rank"]][1])
  } else {
    return(smooth_term[["bs.dim"]] - smooth_term[["null.space.dim"]])
  }
}

generate_smooth_terms <- function(covariates, data) {
  if (length(covariates) == 1) {
    # Univariate smooth
    formula_str <- paste0("s(", covariates[1], ")")
  } else {
    # Bivariate smooth
    formula_str <- paste0("s(", paste(covariates, collapse=", "), ")")
  }
  sm_formula <- eval(parse(text = formula_str))
  sm <- mgcv::smoothCon(sm_formula, absorb.cons=TRUE, data = data)[[1]]
  re <- mgcv::smooth2random(sm, "", type = 2)
  pred_matrix <- s2rPred(sm, re, data = data)
  list(Xf = pred_matrix$Xf, Xr = pred_matrix$rand[[1]], n_knots = calculate_n_knots(sm), type = ifelse(length(covariates) == 1, "univariate", "bivariate"))
}

results <- future_map(covariates, ~generate_smooth_terms(.x, data = temp_data)) 

for (i in seq_along(covariates)) {
  result_type <- results[[i]]$type
  if (result_type == "univariate") {
    covariate <- covariates[[i]]
    temp_data[[paste0("Xf_", covariate)]] <- results[[i]]$Xf
    temp_data[[paste0("Xr_", covariate)]] <- results[[i]]$Xr
    temp_data[[paste0("n_knots_", covariate)]] <- results[[i]]$n_knots
  } else { # bivariate case
    covariate1 <- covariates[[i]][1]
    covariate2 <- covariates[[i]][2]
    bivariate_name <- paste0(covariate1, "_", covariate2)
    temp_data[[paste0("Xf_", bivariate_name)]] <- results[[i]]$Xf
    temp_data[[paste0("Xr_", bivariate_name)]] <- results[[i]]$Xr
    temp_data[[paste0("n_knots_", bivariate_name)]] <- results[[i]]$n_knots
  }
}

n_knots_list <- setNames(lapply(results, `[[`, "n_knots"), 
                         sapply(covariates, function(x) {
                           if (is.character(x) && length(x) == 1) {
                             return(x)
                           } else {
                             return(paste(x, collapse = "_"))
                           }
                         }))

augment_response <- function(y, augmented_length) {
  c(y, rep(0, augmented_length - length(y)))
}

augment_data <- function(data, lambda, n_knots_list) {
  augmented_matrices <- list()
  
  # Calculate the maximum size required for augmentation
  max_augmentation_size <- max(sapply(n_knots_list, function(knots) knots)) + length(lambda)
  
  for (name in names(n_knots_list)) {
    xr_var_name <- paste0("Xr_", name)
    xr_var <- data[[xr_var_name]]
    if (is.null(xr_var)) {
      stop(paste0("Variable ", xr_var_name, " is NULL"))
    }
    
    # Augment each matrix to the maximum size
    augmented_matrix <- augment_design_matrix(xr_var, lambda, max_augmentation_size)
    augmented_matrices[[name]] <- augmented_matrix
  }
  
  # Augment the response variable as well
  y_augmented <- augment_response(data$degrees_from_mean, max_augmentation_size)
  
  augmented_data <- data.frame(degrees_from_mean = y_augmented)
  for (name in names(n_knots_list)) {
    augmented_matrix <- augmented_matrices[[name]]
    augmented_data[[paste0("Xr_", name)]] <- I(augmented_matrix)
  }
  
  return(augmented_data)
}


augment_design_matrix <- function(Xr, lambda, max_augmentation_size) {
  if (is.null(Xr) || nrow(Xr) == 0 || ncol(Xr) == 0) {
    stop("Input matrix Xr is NULL or empty")
  }
  
  # Calculate the number of rows to be added to match the max_augmentation_size
  additional_rows <- max_augmentation_size - nrow(Xr)
  
  # Augment the matrix with additional rows for ridge penalty
  tryCatch({
    augmented_matrix <- rbind(Xr, sqrt(lambda) * diag(ncol(Xr)), matrix(0, nrow = additional_rows, ncol = ncol(Xr)))
    if (is.null(augmented_matrix) || nrow(augmented_matrix) == 0) {
      stop("Augmented matrix is NULL or empty after applying ridge penalty")
    }
  }, error = function(e) {
    stop(paste("Error in augment_design_matrix: ", e$message))
  })
  
  return(augmented_matrix)
}


gcv_lambda <- function(lambda, n_knots_list, data) {
  augmented_data <- augment_data(data, lambda, n_knots_list)
  
  gcv_model <- glmmTMB(formula = degrees_from_mean ~ 1 + Xr_time + Xr_min_temp + Xr_max_temp +
                         + Xr_latitude_longitude, 
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

lambda_values <- seq(0.001, 0.2, by = 0.02)
cv_results <- cross_validate_lambda(temp_data, lambda_values, n_knots_list)
optimal_lambda <- cv_results$optimal_lambda
gcv_values <- cv_results$gcv_scores
temp_data_augmented <- augment_data(temp_data, optimal_lambda, n_knots_list)

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
total_rows <- nrow(temp_data)
train_size <- floor(train_proportion * total_rows)
valid_size <- floor(valid_proportion * total_rows)

train_data <- temp_data[1:train_size, ]
valid_data <- temp_data[(train_size + 1):(train_size + valid_size), ]
test_data <- temp_data[(train_size + valid_size + 1):total_rows, ]

train_data_augmented <- temp_data_augmented[1:train_size, ]
valid_data_augmented <- temp_data_augmented[(train_size + 1):(train_size + valid_size), ]
test_data_augmented <- temp_data_augmented[(train_size + valid_size + 1):total_rows, ]

calculate_metrics <- function(actual, predicted) {
  sqrt(mean((predicted - actual)^2))
}

fit_and_evaluate_model <- function(model_formula, train_data, valid_data, family, use_REML = FALSE) {
  start_time <- Sys.time()
  model <- glmmTMB(model_formula, data = train_data, family = family, REML = use_REML)
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
      AutoModel = fit_and_evaluate_model(auto_model_formula, cv_train, cv_valid, auto_model_family, use_REML = FALSE)
    )
  })
  return(metrics_list)
}

k <- 5
folds <- cut(seq(1, nrow(temp_data)), breaks = k, labels = FALSE)
results <- cross_validate_models(temp_data, temp_data_augmented, folds)
avg_metrics <- map(results, ~sapply(.x, function(x) sapply(x, mean)))
print(avg_metrics)

fit_and_predict <- function(model_formula, model_family, train_data, test_data, augmented = FALSE, augmented_train_data = NULL, use_REML = FALSE) {
  model_data <- if (augmented && !is.null(augmented_train_data)) augmented_train_data else train_data
  model <- glmmTMB(model_formula, data = model_data, family = model_family, REML = use_REML)
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
  ManualModel = list(formula = manual_model_formula, family = manual_model_family, augmented = TRUE, use_REML = FALSE),
  AutoModel = list(formula = auto_model_formula, family = auto_model_family, augmented = FALSE, use_REML = FALSE)
)

model_predictions <- future_map(models_to_train, function(model_info) {
  fit_and_predict(model_info$formula, model_info$family, full_train_data, test_data, augmented = model_info$augmented, augmented_train_data = full_train_data_augmented, use_REML = model_info$use_REML)
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
    geom_smooth(method = "loess", se = TRUE, linewidth = 1, aes(color = Model)) +  # Color matched to Model, fully opaque
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
      color = guide_legend(override.aes = list(alpha = 1, size = 3)),
      shape = guide_legend(override.aes = list(size = 3))
    )
}

plot_predictions(model_predictions, test_data, models_to_train)

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
library(DHARMa)
library(reshape2)
library(purrr)
library(furrr)
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

temp_data <- read_csv("weather-anomalies-1964-2013.csv")

temp_data$date_str <- as.Date(temp_data$date_str)

temp_data <- temp_data %>%
  arrange(date_str) %>%
  mutate(time = as.numeric(date_str - min(date_str)) + 1)

temp_data <- temp_data %>%
  mutate(month = factor(month(date_str), levels = 1:12))


# Filter for years 1990 to 2010 and for locations within the contiguous US bounds
temp_data <- temp_data %>%
  filter(year(date_str) >= 1990 & year(date_str) <= 2010,
         latitude >= 24.396308 & latitude <= 49.384358,
         longitude >= -125.001651 & longitude <= -66.93457)

# Calculate grid sizes and unique stations as before
grid_size_lat <- 2
grid_size_lon <- 2

unique_stations <- temp_data %>%
  mutate(lat_grid = floor(latitude / grid_size_lat),
         lon_grid = floor(longitude / grid_size_lon)) %>%
  group_by(lat_grid, lon_grid) %>%
  summarise(id = first(id)) %>%
  ungroup()

unique_station_ids <- unique(unique_stations$id)

filtered_data <- temp_data %>%
  filter(id %in% unique_station_ids)

temp_data <- filtered_data
rm(filtered_data)

temp_data$hot_cold_indicator <- ifelse(temp_data$degrees_from_mean >= 0, 1, 0)

temp_data$date_str <- NULL
total_rows <- nrow(temp_data)
train_size <- floor(0.6 * total_rows)
valid_size <- floor(0.15 * total_rows)

train_data <- temp_data[1:train_size, ]
valid_data <- temp_data[(train_size + 1):(train_size + valid_size), ]
test_data <- temp_data[(train_size + valid_size + 1):total_rows, ]

n_folds <- 5
min_initial_obs <- floor(0.7 * train_size)
remaining_obs <- train_size - min_initial_obs
adjusted_fold_size <- floor(remaining_obs / (n_folds - 1))

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
  
  # Initialize final_valid_index at the start to ensure it's always defined
  final_valid_index <- min(train_size, current_train_end + adjusted_fold_size)

  if (i < n_folds) {
    current_valid_data <- train_data[(current_train_end + 1):(current_train_end + adjusted_fold_size), ]
  } else {
    current_valid_data <- train_data[(current_train_end + 1):final_valid_index, ]
  }
  
  # Spatial Feature
  coordinates <- current_train_data[, c("longitude", "latitude")]
  coordinates_sp <- SpatialPoints(coords = as.matrix(coordinates))
  nb <- knn2nb(knearneigh(coordinates_sp, k = 5)) # Adjust k as needed
  weights <- nb2listw(nb, style = "W")
  current_train_data$spatial <- lag.listw(weights, current_train_data$degrees_from_mean)
  spatial_forecast <- forecast(current_train_data$spatial, h = nrow(current_valid_data))
  
  ts_data <- ts(current_train_data$degrees_from_mean, frequency = 365) # Adjust frequency as needed
  decomposed <- stl(ts_data, s.window = "periodic", robust = TRUE)
  current_train_data$seasonal <- decomposed$time.series[, "seasonal"]
  current_train_data$trend <- decomposed$time.series[, "trend"]
  
  # Forecast for the length of the current validation set
  trend_forecast <- forecast(current_train_data$trend, h = nrow(current_valid_data))
  seasonal_forecast <- forecast(current_train_data$seasonal, h = nrow(current_valid_data))
  # Add forecasted trend values to current_train_data and current_valid_data
  
  current_valid_data$spatial <- spatial_forecast$mean
  current_valid_data$trend <- trend_forecast$mean  
  current_valid_data$seasonal <- seasonal_forecast$mean
  
  cv_model <- glmmTMB(hot_cold_indicator ~ s(time) + 
                        s(max_temp)*s(min_temp) +
                        s(latitude,longitude) + s(spatial) +
                        s(trend) + s(seasonal) +
                        (1|month), 
                      data = current_train_data,
                      family = binomial(link="logit"),
                      REML = TRUE)
  
  predictions <- predict(cv_model, newdata = current_valid_data, type = "response")
  predicted_classes <- ifelse(predictions >= 0.5, 1, 0)
  
  conf_matrix <- confusionMatrix(as.factor(predicted_classes), as.factor(current_valid_data$hot_cold_indicator))
  f1_score <- conf_matrix[["byClass"]][["F1"]]
  
  total_f1_score <- total_f1_score + f1_score
  fold_count <- fold_count + 1
}

average_f1_score <- total_f1_score / fold_count
print(paste("Average F1 Score:", average_f1_score))

# Combine train and validation sets to form the full training dataset
full_train_data <- rbind(train_data, valid_data)

# Spatial Feature
coordinates <- full_train_data[, c("longitude", "latitude")]
coordinates_sp <- SpatialPoints(coords = as.matrix(coordinates))
nb <- knn2nb(knearneigh(coordinates_sp, k = 5)) # Adjust k as needed
weights <- nb2listw(nb, style = "W")
full_train_data$spatial <- lag.listw(weights, full_train_data$degrees_from_mean)

ts_data <- ts(full_train_data$degrees_from_mean, frequency = 5) # Adjust frequency as needed
decomposed <- stl(ts_data, s.window = "periodic", robust = TRUE)
full_train_data$seasonal <- decomposed$time.series[, "seasonal"]
full_train_data$trend <- decomposed$time.series[, "trend"]

spatial_forecast <- forecast(full_trian_data$spatial, h = nrow(test_data))
trend_forecast <- forecast(full_train_data$trend, h = nrow(test_data))
seasonal_forecast <- forecast(full_train_data$seasonal, h = nrow(test_data))

test_data$spatial <- spatial_forecast$mean
test_data$trend <- trend_forecast$mean  
test_data$seasonal <- seasonal_forecast$mean

# Re-train the model on the full training dataset
final_model <- glmmTMB(hot_cold_indicator ~ s(time, k = 5) + s(max_temp, k = 6)*s(min_temp, k = 6) +
                       s(latitude,longitude) + s(trend) + s(spatial) + s(seasonal)
                       + (1|month), 
                       data = full_train_data,
                       family = binomial(link="logit"),
                       REML = TRUE)

# Predict on the unseen test data
test_predictions <- predict(final_model, newdata = test_data, type = "response")
test_predicted_classes <- ifelse(test_predictions >= 0.5, 1, 0)


test_conf_matrix <- confusionMatrix(as.factor(test_predicted_classes), as.factor(test_data$hot_cold_indicator))
test_f1_score <- test_conf_matrix[["byClass"]][["F1"]]
print(paste("Test F1 Score:", test_f1_score))

F1_score <- test_conf_matrix[["byClass"]][["F1"]]
print(paste("F1 Score:", F1_score))

conf_matrix_table <- test_conf_matrix$table
conf_matrix_df <- as.data.frame(conf_matrix_table)
names(conf_matrix_df) <- c("Reference", "Prediction", "Freq")

ggplot(data = conf_matrix_df, aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%d", Freq)), vjust = 1) +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal() +
  labs(title = "Confusion Matrix", x = "Actual Class", y = "Predicted Class")



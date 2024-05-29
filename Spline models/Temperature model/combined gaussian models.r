source("requirements.R")

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
# Base map of the USA 
usa_map <- map_data("state")
usa_map <- usa_map[!usa_map$region %in% c("alaska", "hawaii"),]
# Check the dist of our stations
ggplot() +
  geom_polygon(data = usa_map, aes(x = long, y = lat, group = group), 
               fill = "lightgrey", color = "black") +
  geom_point(data = temp_data, aes(x = longitude, y = latitude), 
             color = "red", size = 1, alpha = 0.5) +
  coord_fixed(1.3) + # This helps to keep the aspect ratio of the USA map
  labs(title = "Geographical Distribution of Weather Stations",
       x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

rm(usa_map)
gc()

temp_data$hot_cold_indicator <- ifelse(temp_data$degrees_from_mean >= 0, 1, 0)

temp_data$date_str <- NULL
temp_data <- temp_data[order(temp_data$time), ]
total_rows <- nrow(temp_data)
train_size <- floor(0.6 * total_rows) 
valid_size <- floor(0.2 * total_rows) 
test_size <- total_rows - train_size - valid_size 

train_data <- temp_data[1:train_size, ]
valid_data <- temp_data[(train_size + 1):(train_size + valid_size), ]
test_data <- temp_data[(train_size + valid_size + 1):total_rows, ]


# Further split based on temperature
train_hot <- train_data[train_data$degrees_from_mean > 0, ]
train_cold <- train_data[train_data$degrees_from_mean <= 0, ]
valid_hot <- valid_data[valid_data$degrees_from_mean > 0, ]
valid_cold <- valid_data[valid_data$degrees_from_mean <= 0, ]
test_hot <- test_data[test_data$degrees_from_mean > 0, ]
test_cold <- test_data[test_data$degrees_from_mean <= 0, ]

calculate_rmse <- function(actual, predicted) {
  sqrt(mean((actual - predicted)^2))
}

n_folds <- 5
min_initial_obs <- floor(0.8 * train_size)
remaining_obs <- train_size - min_initial_obs
adjusted_fold_size <- floor(remaining_obs / (n_folds - 1))

hot_rmse_values <- c()
cold_rmse_values <- c()

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
    current_valid_data <- train_data[(current_train_end + 1):train_size, ]
  }
  
  # Spatial Feature
  coordinates <- current_train_data[, c("longitude", "latitude")]
  coordinates_sp <- SpatialPoints(coords = as.matrix(coordinates))
  nb <- knn2nb(knearneigh(coordinates_sp, k = 5)) # Adjust k as needed
  weights <- nb2listw(nb, style = "W")
  current_train_data$spatial <- lag.listw(weights, current_train_data$degrees_from_mean)
  spatial_forecast <- forecast(current_train_data$spatial, h = nrow(current_valid_data))
  
  ts_data <- ts(current_train_data$degrees_from_mean, frequency = 5) # Adjust frequency as needed
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
  
  current_train_hot <- current_train_data[current_train_data$degrees_from_mean > 0, ]
  current_train_cold <- current_train_data[current_train_data$degrees_from_mean <= 0, ]

  hotmod_cv <- glmmTMB(degrees_from_mean ~ s(time) + s(max_temp)*s(min_temp) +
                            s(latitude)*s(longitude) + 
                            s(trend) + s(spatial) + s(seasonal) +
                            (1|month) + (1|type), 
                         disp =~ 1,
                         family = gaussian(link="identity"), 
                         data = current_train_hot, 
                         REML = TRUE)
  
  coldmod_cv <- glmmTMB(degrees_from_mean ~ s(time) + s(max_temp)*s(min_temp) +
                             s(latitude)*s(longitude) + 
                             s(trend) + s(spatial) + s(seasonal) +
                             (1|month) + (1|type),
                          disp =~ 1,
                          family = gaussian(link="identity"), 
                          data = current_train_cold, 
                          REML = TRUE)
  
  current_valid_hot <- current_valid_data[current_valid_data$degrees_from_mean > 0, ]
  current_valid_cold <- current_valid_data[current_valid_data$degrees_from_mean <= 0, ]
  
  hot_predictions <- predict(hotmod_cv, current_valid_hot, type="response", allow.new.levels = TRUE)
  cold_predictions <- predict(coldmod_cv, current_valid_cold, type="response", allow.new.levels = TRUE)

  hot_rmse <- calculate_rmse(current_valid_hot$degrees_from_mean, hot_predictions)
  cold_rmse <- calculate_rmse(current_valid_cold$degrees_from_mean, cold_predictions)
  
  hot_rmse_values <- c(hot_rmse_values, hot_rmse)
  cold_rmse_values <- c(cold_rmse_values, cold_rmse)
  
  cat(sprintf("Fold %d: Hot Model RMSE = %f, Cold Model RMSE = %f\n", i, hot_rmse, cold_rmse))
  
  pb$tick()
}

average_hot_rmse <- mean(hot_rmse_values)
average_cold_rmse <- mean(cold_rmse_values)

cat(sprintf("Average Hot Model RMSE: %f\n", average_hot_rmse))
cat(sprintf("Average Cold Model RMSE: %f\n", average_cold_rmse))


train_hot_full <- rbind(train_hot, valid_hot)
train_cold_full <- rbind(train_cold, valid_cold)

# Spatial Feature
coordinates <- full_train_hot[, c("longitude", "latitude")]
coordinates_sp <- SpatialPoints(coords = as.matrix(coordinates))
nb <- knn2nb(knearneigh(coordinates_sp, k = 5)) # Adjust k as needed
weights <- nb2listw(nb, style = "W")
full_train_hot$spatial <- lag.listw(weights, full_train_hot$degrees_from_mean)

ts_data <- ts(full_train_hot$degrees_from_mean, frequency = 5) # Adjust frequency as needed
decomposed <- stl(ts_data, s.window = "periodic", robust = TRUE)
full_train_hot$seasonal <- decomposed$time.series[, "seasonal"]
full_train_hot$trend <- decomposed$time.series[, "trend"]

spatial_forecast <- forecast(full_trian_hot$spatial, h = nrow(test_hot))
trend_forecast <- forecast(full_train_hot$trend, h = nrow(test_hot))
seasonal_forecast <- forecast(full_train_hot$seasonal, h = nrow(test_hot))

test_hot$spatial <- spatial_forecast$mean
test_hot$trend <- trend_forecast$mean  
test_hot$seasonal <- seasonal_forecast$mean

# Spatial Feature
coordinates <- full_train_cold[, c("longitude", "latitude")]
coordinates_sp <- SpatialPoints(coords = as.matrix(coordinates))
nb <- knn2nb(knearneigh(coordinates_sp, k = 5)) # Adjust k as needed
weights <- nb2listw(nb, style = "W")
full_train_cold$spatial <- lag.listw(weights, full_train_cold$degrees_from_mean)

ts_data <- ts(full_train_cold$degrees_from_mean, frequency = 5) # Adjust frequency as needed
decomposed <- stl(ts_data, s.window = "periodic", robust = TRUE)
full_train_cold$seasonal <- decomposed$time.series[, "seasonal"]
full_train_cold$trend <- decomposed$time.series[, "trend"]

spatial_forecast <- forecast(full_trian_cold$spatial, h = nrow(test_cold))
trend_forecast <- forecast(full_train_cold$trend, h = nrow(test_cold))
seasonal_forecast <- forecast(full_train_cold$seasonal, h = nrow(test_cold))

test_cold$spatial <- spatial_forecast$mean
test_cold$trend <- trend_forecast$mean  
test_cold$seasonal <- seasonal_forecast$mean

hotmod <- glmmTMB(degrees_from_mean ~ s(time) + s(max_temp)*s(min_temp) +
                    s(latitude)*s(longitude) + 
                    s(trend) + s(spatial) + s(seasonal), 
                  disp =~ type,
                  family = gaussian(link="identity"),
                  data = train_hot_full, 
                  REML = TRUE)

coldmod <- glmmTMB(degrees_from_mean ~ s(time) + s(max_temp)*s(min_temp) +
                     s(latitude)*s(longitude) + 
                     s(trend) + s(spatial) + s(seasonal),
                   disp =~ type,
                   family = gaussian(link="identity"),
                   data = train_cold_full, 
                   REML = TRUE)

predictions_hot <- predict(hotmod, newdata = test_hot, type = "response", allow.new.levels=TRUE)
rmse_hot <- sqrt(mean((test_hot$degrees_from_mean - predictions_hot)^2))
predictions_cold <- predict(coldmod, newdata = test_cold, type = "response", allow.new.levels=TRUE)
rmse_cold <- sqrt(mean((test_cold$degrees_from_mean - predictions_cold)^2))
print(paste("RMSE for hot deviations:", rmse_hot))
print(paste("RMSE for cold deviations:", rmse_cold))

test_hot$predicted_degrees_from_mean <- predictions_hot
test_cold$predicted_degrees_from_mean <- predictions_cold

test_data$predicted_degrees_from_mean[test_data$degrees_from_mean > 0] <- test_hot$predicted_degrees_from_mean
test_data$predicted_degrees_from_mean[test_data$degrees_from_mean <= 0] <- test_cold$predicted_degrees_from_mean
test_data$actual <- test_data$degrees_from_mean
test_data$predicted <- test_data$predicted_degrees_from_mean
rmse_combined <- sqrt(mean((test_data$actual - test_data$predicted)^2, na.rm = TRUE))
print(rmse_combined)

ggplot(test_data, aes(x = actual, y = predicted)) +
  geom_point(alpha = 0.3, color = "lightblue") + 
  geom_smooth(method = "loess", se = FALSE, color = "darkblue") +   
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") + 
  labs(title = "Predicted vs Actual Degrees From Mean", 
       x = "Actual Degrees From Mean", 
       y = "Predicted Degrees From Mean") +
  theme_minimal() 

long_data <- reshape2::melt(test_data, id.vars = "id", measure.vars = c("predicted", "actual"))

ggplot(long_data, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.6) +  
  labs(x = "Degrees From Mean", y = "Density", 
       title = "Density Plot of Predicted vs. Actual Degrees From Mean") +
  scale_fill_manual(values = c("predicted" = "darkblue", "actual" = "lightblue"),
                    labels = c("Predicted", "Actual")) +
  theme_minimal()



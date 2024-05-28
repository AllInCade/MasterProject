library(readr)
library(readxl)
library(dplyr)
library(lubridate)
library(mgcv)
library(glmmTMB)
library(gamm4)
library(forecast)
library(caret)
library(ggplot2)
library(DHARMa)
library(reshape2)
library(MASS)
library(car)
library(VGAM)
library(optimx)
library(bestNormalize)
library(randomForest)
# DATA
# https://banks.data.fdic.gov/explore/failures?aggReport=detail&displayFields=NAME%2CCERT%2CFIN%2CCITYST%2CFAILDATE%2CSAVR%2CRESTYPE%2CCOST%2CRESTYPE1%2CCHCLASS1%2CQBFDEP%2CQBFASSET&endFailYear=2024&sortField=CERT&sortOrder=asc&startFailYear=1975

bank_data <- read_csv("bankdata.csv")
bank_data <- bank_data[bank_data$ID != 613, ]
bank_data <- na.omit(bank_data)
bank_data$FAILDATE <- as.Date(bank_data$FAILDATE, format = "%m/%d/%Y")

split_city_state <- strsplit(bank_data$CITYST, ", ")
bank_data$CITY <- sapply(split_city_state, function(x) x[1])
bank_data$STATE <- sapply(split_city_state, function(x) ifelse(length(x) > 1, x[2], NA))

bank_data <- bank_data %>%
  mutate(
    YEAR = year(FAILDATE),
    MONTH = month(FAILDATE),
    WEEKDAY = wday(FAILDATE, label = TRUE)
  )

start_date <- min(bank_data$FAILDATE, na.rm = TRUE)
bank_data$TIME <- as.integer(difftime(bank_data$FAILDATE, start_date, units = "days")) + 1
bank_data$QBFASSET_scaled <- (bank_data$QBFASSET - median(bank_data$QBFASSET, na.rm = TRUE)) / IQR(bank_data$QBFASSET, na.rm = TRUE)
bank_data$COST_scaled <- (bank_data$COST - median(bank_data$COST, na.rm = TRUE)) / IQR(bank_data$COST, na.rm = TRUE)
bank_data$DEP_ASS_RATIO <- bank_data$QBFDEP + 1e-6/ (bank_data$QBFASSET)
bank_data$DEP_ASS_RATIO[is.infinite(bank_data$DEP_ASS_RATIO)] <- NA
bank_data$DEP_ASS_RATIO[is.nan(bank_data$DEP_ASS_RATIO)] <- NA
# Predict COST 
# Feature ImportanceAnalysis
bank_data <- bank_data[bank_data$COST > 0, ]

## Due to extreme outliers and the extremely skewness and kurtosis of the data 
## we will log transform numerical predictors to stabilize 

bank_data$log_COST <- log(bank_data$COST + 1)
bank_data$log_ASSET <- log(bank_data$QBFASSET + 1)
bank_data$log_DEPOSIT <- log(bank_data$QBFDEP + 1)

write_csv(bank_data, "bank_data_r.csv")

set.seed(101) 
training_indices <- sample(1:nrow(bank_data), 0.8 * nrow(bank_data))
train_data <- bank_data[training_indices, ]
test_data <- bank_data[-training_indices, ]

model <- randomForest(log_COST ~ . - COST - COST_scaled - log_COST - CITYST - FAILDATE - CERT - FIN - ID - QBFASSET_scaled - QBFASSET - QBFDEP - DEP_ASS_RATIO, data = train_data, na.action = na.omit, importance = TRUE)

predictions <- predict(model, test_data)
actuals <- test_data$log_COST

RMSE <- sqrt(mean((predictions - actuals)^2))
print(paste("RMSE:", RMSE))



importance_scores <- importance(model, type = 1) 
feature_names <- row.names(importance_scores)
importance_data <- data.frame(Feature = feature_names, Importance = importance_scores[, 1])
sorted_importance_data <- importance_data[order(-importance_data$Importance), ] # Sort in decreasing order
barplot(sorted_importance_data$Importance, names.arg = sorted_importance_data$Feature, 
        main = "Feature Importance in Predicting COST", las = 2, cex.names = 0.7)


## check if smooths can offer increased performance

linmod <- glmmTMB(log_COST ~ TIME + log_DEPOSIT +log_ASSET, 
                  data = bank_data)
  
tmbmod <- glmmTMB(log_COST ~ s(TIME) + s(log_ASSET) + s(log_DEPOSIT),
                  data = bank_data, 
                  family = gaussian(link = "identity"),
                  REML = TRUE)

AIC(linmod)
AIC(tmbmod) 

# GAM shows better in-sample (as one may expect) performance in termf of AIC
# Explore more advanced options like random effects etc

tmbmod2 <- glmmTMB(log_COST ~ s(TIME) + s(log_ASSET) + s(log_DEPOSIT)
                   + (1|YEAR) + (1|STATE) + (1|CHCLASS1),
                   dispformula =~ SAVR + RESTYPE,
                   data = bank_data, 
                   family = gaussian(),
                   REML = TRUE)



gammod2 <- gamm4(log_COST ~ s(TIME) + s(log_ASSET) + s(log_DEPOSIT),
                   random =~ (1|YEAR) + (1|STATE) + (1|CHCLASS1),
                   data = bank_data, 
                   family = gaussian(link = "identity"), 
                   REML = TRUE)



AIC(tmbmod2)
AIC(gammod2$mer)

set.seed(420)  

n <- nrow(bank_data)
train_size <- round(0.6 * n)
valid_size <- round(0.2 * n)
bank_data_shuffled <- bank_data[sample(nrow(bank_data)), ]
train_data <- bank_data_shuffled[1:train_size, ]
valid_data <- bank_data_shuffled[(train_size + 1):(train_size + valid_size), ]
test_data <- bank_data_shuffled[(train_size + valid_size + 1):n, ]


tmbmod2_train <- glmmTMB(log_COST ~ s(TIME) + s(log_ASSET) + s(log_DEPOSIT)
                         + (1|YEAR) + (1|STATE) + (1|CHCLASS1),
                         dispformula =~ SAVR + RESTYPE,
                         data = train_data, 
                         family = gaussian(),
                         REML = TRUE)

gammod2_train <-gamm4(log_COST ~ s(TIME) + s(log_ASSET) + s(log_DEPOSIT),
                        random =~ (1|YEAR) + (1|STATE) + (1|CHCLASS1),
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

#valid_data$WEEKDAY <- factor(valid_data$WEEKDAY, levels = sort(unique(valid_data$WEEKDAY)))
# Loop through each fold
for(i in 1:n_folds) {
  # Define training and test sets for this fold
  test_indices <- which(folds == i)
  train_indices <- setdiff(seq_len(nrow(valid_data)), test_indices)
  train_fold <- valid_data[train_indices, ]
  test_fold <- valid_data[test_indices, ]
  
  model_gamm <- gamm4(log_COST ~ s(TIME) + s(log_ASSET) + s(log_DEPOSIT),
                      random =~ (1|YEAR) + (1|STATE) + (1|CHCLASS1),
                      data = train_fold, 
                      family = gaussian(link = "identity"), 
                      REML = TRUE)
  pred_gamm <- predict(model_gamm$gam, test_fold, type = "response")
  rmse_gamm[[i]] <- sqrt(mean((pred_gamm - test_fold$log_COST)^2))
  
  model_glmmTMB <- glmmTMB(log_COST ~ s(TIME) + s(log_ASSET) + s(log_DEPOSIT)
                           + (1|YEAR) + (1|STATE) + (1|CHCLASS1) + (1|RESTYPE1),
                           dispformula =~ RESTYPE,
                           data = train_fold, 
                           family = t_family("identity"),
                           start = list(psi = log(4)),
                           REML = TRUE)
  pred_glmmTMB <- predict(model_glmmTMB, test_fold, type = "response", allow.new.levels=TRUE)
  rmse_glmmTMB[[i]] <- sqrt(mean((pred_glmmTMB - test_fold$log_COST)^2))
}

mean_rmse_gamm <- mean(unlist(rmse_gamm), na.rm = TRUE) # Using na.rm=TRUE to handle potential NA values
mean_rmse_glmmTMB <- mean(unlist(rmse_glmmTMB), na.rm = TRUE)
cat("Average RMSE for GAMM: ", mean_rmse_gamm, "\n")
cat("Average RMSE for GLMMTMB: ", mean_rmse_glmmTMB, "\n")


full_train_data <- rbind(train_data, valid_data)
start_time_gamm_train <- Sys.time()
full_model_gamm <- model_gamm <- gamm4(log_COST ~ s(TIME) + s(log_ASSET) + s(log_DEPOSIT),
                                       random =~ (1|YEAR) + (1|STATE) + (1|CHCLASS1) + (1|RESTYPE1),
                                       data = full_train_data, 
                                       family = gaussian(link = "identity"), 
                                       REML = TRUE)

end_time_gamm_train <- Sys.time()
train_duration_gamm <- end_time_gamm_train - start_time_gamm_train
start_time_gamm_pred <- Sys.time()
predictions_gamm <- predict(full_model_gamm$gam, test_data, type = "response")
end_time_gamm_pred <- Sys.time()
pred_duration_gamm <- end_time_gamm_pred - start_time_gamm_pred
rmse_gamm <- sqrt(mean((predictions_gamm - test_data$log_COST)^2))

start_time_glmmTMB_train <- Sys.time()

full_model_glmmTMB <- glmmTMB(log_COST ~ s(TIME) + s(log_ASSET) + s(log_DEPOSIT)
                              + (1|YEAR) + (1|STATE) + (1|CHCLASS1) + (1|RESTYPE1),
                              dispformula =~ SAVR + RESTYPE,
                              data = full_train_data, 
                              family = t_family("identity"),
                              start = list(psi = log(4)),
                              REML = TRUE)

  end_time_glmmTMB_train <- Sys.time()
train_duration_glmmTMB <- end_time_glmmTMB_train - start_time_glmmTMB_train

start_time_glmmTMB_pred <- Sys.time()
predictions_glmmTMB <- predict(full_model_glmmTMB, test_data, type = "response", allow.new.levels=TRUE)
end_time_glmmTMB_pred <- Sys.time()
pred_duration_glmmTMB <- end_time_glmmTMB_pred - start_time_glmmTMB_pred
rmse_glmmTMB <- sqrt(mean((predictions_glmmTMB - test_data$log_COST)^2))

cat("Training duration for GAMM: ", train_duration_gamm, "\n")
cat("Training duration for GLMMTMB: ", train_duration_glmmTMB, "\n")

cat("Prediction duration for GAMM: ", pred_duration_gamm, "\n")
cat("Prediction duration for GLMMTMB: ", pred_duration_glmmTMB, "\n")

# Print the results
cat("RMSE for GAMM: ", rmse_gamm, "\n")
cat("RMSE for GLMMTMB: ", rmse_glmmTMB, "\n")

# Correctly combine predictions for the plot
plot_data <- data.frame(
  Actual = rep(test_data$log_COST, 2),
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


actual_density <- density(full_train_data$log_COST)

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
  scale_color_manual(values = c("Actual" = "black")) + # Assign color to the actual values density line
  scale_fill_manual(values = c("GAMM" = "darkorange", "GLMMTMB" = "darkblue")) +
  labs(title = "Density of Predicted Values vs Actual Values",
       x = "Values",
       y = "Density",
       fill = "Model",
       color = "Density Type") + # Updated to include 'Actual' label
  theme_minimal() +
  theme(legend.position = "right") +
  facet_wrap(~Model)

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
#install.packages("clock")
#install.packages("recipes")
#install.packages("bestNormalize")
library(bestNormalize)

# DATA
# https://banks.data.fdic.gov/explore/failures?aggReport=detail&displayFields=NAME%2CCERT%2CFIN%2CCITYST%2CFAILDATE%2CSAVR%2CRESTYPE%2CCOST%2CRESTYPE1%2CCHCLASS1%2CQBFDEP%2CQBFASSET&endFailYear=2024&sortField=CERT&sortOrder=asc&startFailYear=1975

bank_data_full <- read_csv("bankdata.csv")
bank_data_full <- bank_data_full[bank_data_full$ID != 613, ]
bank_data_na_omit <- na.omit(bank_data_full)
bank_data_na_omit$FAILDATE <- as.Date(bank_data_na_omit$FAILDATE, format = "%m/%d/%Y")

# We'll train a random Forest on the data set with NAs omitted using the assumed
# informative predictors of total assets and total deposits to predict the 
# COST variable (Estimated Loss associated with the bank failure).
# Use the model to impute the missing values in the full data set

train_data <- bank_data_na_omit[, c("COST", "QBFASSET", "QBFDEP")]

control <- trainControl(method="cv", number=10)
tuneGrid <- expand.grid(.mtry=c(1, 2, round(sqrt(ncol(train_data)-1)), ncol(train_data)-1))

set.seed(123) 
model_cv <- train(COST ~ ., data=train_data, method="rf",
                  trControl=control, tuneGrid=tuneGrid,
                  importance=TRUE)

print(model_cv)

if(!is.null(model_cv$bestTune)) {
  best_mtry <- model_cv$bestTune$mtry
} else {
  stop("No valid 'bestTune' object found.")
}

if(!is.integer(best_mtry) || best_mtry < 1 || best_mtry > ncol(train_data) - 1) {
  stop("The 'mtry' value is invalid or out of the valid range.")
}

if(!is.null(best_mtry)) {
  final_model <- randomForest(COST ~ ., data=train_data, mtry=best_mtry, importance=TRUE)
} else {
  stop("The 'mtry' value is NULL.")
}

ratios <- with(bank_data_na_omit, QBFASSET / QBFDEP)
median_ratio <- median(ratios, na.rm = TRUE)

# Find the index of the row with missing QBFASSET and known QBFDEP
missing_qbfasset_index <- which(is.na(bank_data_full$QBFASSET) & !is.na(bank_data_full$QBFDEP))
bank_data_full$QBFASSET[missing_qbfasset_index] <- median_ratio * bank_data_full$QBFDEP[missing_qbfasset_index]
missing_cost_indices <- which(is.na(bank_data_full$COST))

# Make sure there are no longer any missing values in the predictors
if (any(is.na(bank_data_full[missing_cost_indices, c("QBFASSET", "QBFDEP")]))) {
  stop("predict_data still contains missing values in the predictor columns.")
}

predict_data <- bank_data_full[missing_cost_indices, c("QBFASSET", "QBFDEP")]
predicted_cost <- predict(final_model, newdata = predict_data)
bank_data_full$COST[missing_cost_indices] <- predicted_cost

# Check for NAs in each column
na_count_per_column <- sapply(bank_data_full, function(x) sum(is.na(x)))
total_na_count <- sum(is.na(bank_data_full))
print(na_count_per_column)
print(paste("Total NAs in the dataframe:", total_na_count))

# Prepare data for daily counts of bank failures
bank_data_full$FAILDATE <- as.Date(bank_data_full$FAILDATE, format="%m/%d/%Y")
daily_failures <- aggregate(ID ~ FAILDATE, data = bank_data_full, FUN = length)
colnames(daily_failures)[which(colnames(daily_failures) == "ID")] <- "FailuresCount"
bank_data_full <- merge(bank_data_full, daily_failures, by="FAILDATE", all.x=TRUE)

bank_data_full <- bank_data_full %>%
  mutate(
    YEAR = year(FAILDATE),
    MONTH = month(FAILDATE),
    WEEKDAY = wday(FAILDATE, label = TRUE)
  )


# Split the CITYST column into two columns: CITY and STATE
split_city_state <- strsplit(bank_data_full$CITYST, ", ")
bank_data_full$CITY <- sapply(split_city_state, function(x) x[1])
bank_data_full$STATE <- sapply(split_city_state, function(x) ifelse(length(x) > 1, x[2], NA))
head(bank_data_full[c("CITYST", "CITY", "STATE")])

start_date <- min(bank_data_full$FAILDATE, na.rm = TRUE)
bank_data_full$TIME <- as.integer(difftime(bank_data_full$FAILDATE, start_date, units = "days")) + 1

tmb_pois <- glmmTMB(FailuresCount ~ s(TIME)  + (1|STATE),
                    data = bank_data_full, 
                    family = poisson,
                    REML = TRUE)


gamm_pois <- gamm4(FailuresCount ~ s(TIME), 
                   random = ~(1 | STATE),
                   data = bank_data_full, 
                   family = poisson,
                   REML = TRUE) 


summary(tmb_pois)
summary(gamm_pois$mer)

AIC(tmb_pois)
AIC(gamm_pois$mer)

# The models are approximately equal, but some differences in the optimization
# methods are very difficult to equalize exactly. gamm4 is not easily managed 
# through control options

# Trying neg bin models 
tmb_nb <- glmmTMB(FailuresCount ~ s(TIME) + s(YEAR) + s(MONTH)
                    + (1|STATE) + (1|SAVR) + (1|CHCLASS1) + (1|RESTYPE1) + (1|WEEKDAY),
                    dispformula =~ RESTYPE,
                    data = bank_data_full, 
                    family = nbinom2(),
                    REML = TRUE)

gamm_nb <- gamm4(FailuresCount ~ s(TIME) + s(YEAR) + s(MONTH), 
                   random =~ (1|SAVR) + (1|CHCLASS1) + (1|RESTYPE1),
                   data = bank_data_full, 
                   family = negbin(theta=2),
                   REML = TRUE) 


AIC(tmb_nb)
AIC(gamm_nb$mer)

bank_data <- bank_data_full
# Note data is of type "event-log" and not TS, so random sampling is okay. 
set.seed(420)  

n <- nrow(bank_data)
train_size <- round(0.6 * n)
valid_size <- round(0.2 * n)
bank_data_shuffled <- bank_data[sample(nrow(bank_data)), ]
train_data <- bank_data_shuffled[1:train_size, ]
valid_data <- bank_data_shuffled[(train_size + 1):(train_size + valid_size), ]
test_data <- bank_data_shuffled[(train_size + valid_size + 1):n, ]

tmb_nb_train <- glmmTMB(FailuresCount ~ s(TIME) + s(YEAR) + s(MONTH)
                        + (1|STATE) + (1|SAVR) + (1|CHCLASS1) + (1|RESTYPE1),
                        dispformula =~ RESTYPE + WEEKDAY,
                        data = train_data, 
                        family = nbinom2(link="inverse"),
                        REML = TRUE)

gamm_nb_train <- gamm4(FailuresCount ~ s(TIME) + s(YEAR) + s(MONTH), 
                       random =~ (1|SAVR) + (1|CHCLASS1) + (1|RESTYPE1),
                       data = train_data, 
                       family = negbin(theta=2),
                       REML = TRUE) 

AIC(tmb_nb_train)
AIC(gamm_nb_train$mer)

n_folds <- 5
folds <- cut(seq(1, nrow(valid_data)), breaks = n_folds, labels = FALSE)

rmse_gamm <- vector("list", n_folds)
rmse_glmmTMB <- vector("list", n_folds)

for(i in 1:n_folds) {
  test_indices <- which(folds == i)
  train_indices <- setdiff(seq_len(nrow(valid_data)), test_indices)
  train_fold <- valid_data[train_indices, ]
  test_fold <- valid_data[test_indices, ]
  
  model_gamm <- gamm4(FailuresCount ~ s(TIME) + s(YEAR) + s(MONTH), 
                      random =~ (1|STATE) + (1|SAVR) + (1|CHCLASS1) + (1|WEEKDAY) + (1|RESTYPE),
                      data = train_fold, 
                      family = negbin(theta=2),
                      REML = TRUE)
  
  pred_gamm <- predict(model_gamm$gam, test_fold, type = "response")
  rmse_gamm[[i]] <- sqrt(mean((pred_gamm - test_fold$FailuresCount)^2))
  
  model_glmmTMB <- glmmTMB(FailuresCount ~ s(TIME) + s(YEAR) + s(MONTH)
                           + (1|STATE) + (1|SAVR) + (1|CHCLASS1) + (1|RESTYPE1),
                           dispformula =~ RESTYPE + WEEKDAY,
                           data = train_fold, 
                           family = nbinom2(link="inverse"),
                           REML = TRUE)

  pred_glmmTMB <- predict(model_glmmTMB, test_fold, type = "response", allow.new.levels=TRUE)
  rmse_glmmTMB[[i]] <- sqrt(mean((pred_glmmTMB - test_fold$FailuresCount)^2))
}

mean_rmse_gamm <- mean(unlist(rmse_gamm), na.rm = TRUE) # Using na.rm=TRUE to handle potential NA values
mean_rmse_glmmTMB <- mean(unlist(rmse_glmmTMB), na.rm = TRUE)

cat("Average RMSE for GAMM: ", mean_rmse_gamm, "\n")
cat("Average RMSE for GLMMTMB: ", mean_rmse_glmmTMB, "\n")


full_train_data <- rbind(train_data, valid_data)

start_time_gamm_train <- Sys.time()
full_model_gamm <- gamm4(FailuresCount ~ s(TIME) + s(YEAR) + s(MONTH), 
                         random =~ (1|STATE) + (1|SAVR) + (1|CHCLASS1) + (1|WEEKDAY) + (1|RESTYPE),
                         data = full_train_data, 
                         family = negbin(theta=2),
                         REML = TRUE)

end_time_gamm_train <- Sys.time()
train_duration_gamm <- end_time_gamm_train - start_time_gamm_train
start_time_gamm_pred <- Sys.time()
predictions_gamm <- predict(full_model_gamm$gam, test_data, type = "response")
end_time_gamm_pred <- Sys.time()
pred_duration_gamm <- end_time_gamm_pred - start_time_gamm_pred

fixed_effects <- fixef(tmb_nb_train)
conditional_beta <- fixed_effects[["cond"]]

start_time_glmmTMB_train <- Sys.time()

full_model_glmmTMB <- glmmTMB(FailuresCount ~ s(TIME) + s(YEAR) + s(MONTH)
                              + (1|STATE) + (1|SAVR) + (1|CHCLASS1) + (1|RESTYPE1),
                              dispformula =~ RESTYPE + WEEKDAY,
                              data = full_train_data, 
                              family = nbinom2(link="inverse"),
                              REML = TRUE)

end_time_glmmTMB_train <- Sys.time()
train_duration_glmmTMB <- end_time_glmmTMB_train - start_time_glmmTMB_train

start_time_glmmTMB_pred <- Sys.time()
predictions_glmmTMB <- predict(full_model_glmmTMB, test_data, type = "response", allow.new.levels=TRUE)
end_time_glmmTMB_pred <- Sys.time()
pred_duration_glmmTMB <- end_time_glmmTMB_pred - start_time_glmmTMB_pred

cat("Training duration for GAMM: ", train_duration_gamm, "\n")
cat("Training duration for GLMMTMB: ", train_duration_glmmTMB, "\n")

cat("Prediction duration for GAMM: ", pred_duration_gamm, "\n")
cat("Prediction duration for GLMMTMB: ", pred_duration_glmmTMB, "\n")

plot_data <- data.frame(
  Actual = rep(test_data$FailuresCount, 2),
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
  theme(legend.position = "right") 

actual_density <- density(full_train_data$FailuresCount)

actual_density_df <- data.frame(Value = actual_density$x, Density = actual_density$y)

ggplot() +
  geom_line(data = actual_density_df, aes(x = Value, y = Density, color = "Actual"), size = 1) +
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

rmse_gamm <- sqrt(mean((predictions_gamm - test_data$FailuresCount)^2))
rmse_glmmTMB <- sqrt(mean((predictions_glmmTMB - test_data$FailuresCount)^2))

cat("RMSE for GAMM: ", rmse_gamm, "\n")
cat("RMSE for GLMMTMB: ", rmse_glmmTMB, "\n")

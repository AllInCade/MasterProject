source("requirements.R")


data(ustermlife)
?ustermlife
view(ustermlife)
data=na.omit(ustermlife)
head(data)

summary(data$Face)



# Finding most important variables

set.seed(101) 
training_indices <- sample(1:nrow(data), 0.8 * nrow(data))
train_data <- data[training_indices, ]
test_data <- data[-training_indices, ]

mean(data$Face)

model <- randomForest(Face ~ Age + Gender + Education + MarStat +Ethnicity +SmarStat +Sgender+ Sage +Seducation+ NumHH +Income +TotIncome +Charity +FaceCVLifePol +CashCVLifePol+BorrowCVLifePol +NetValue, data = train_data, na.action = na.omit, importance = TRUE)

predictions <- predict(model, test_data)
actuals <- test_data$Face



RMSE <- sqrt(mean((predictions - actuals)^2))
print(paste("RMSE:", RMSE))



importance_scores <- randomForest::importance(model,type=2) 
feature_names <- row.names(importance_scores)
importance_data <- data.frame(Feature = feature_names, Importance = importance_scores[, 1])
sorted_importance_data <- importance_data[order(-importance_data$Importance), ] # Sort in decreasing order
barplot(sorted_importance_data$Importance, names.arg = sorted_importance_data$Feature, 
        main = "Feature Importance in Predicting Face Value", las = 2, cex.names = 0.7)

print(sorted_importance_data)


insurance_data<-data

insurance_data$Face_log<-log(insurance_data$Face+1)
insurance_data$Income_log<-log(insurance_data$Income+1)
insurance_data$TotIncome_log<-log(insurance_data$TotIncome+1)

insurance_data<-insurance_data[insurance_data$Face < 500000 ,]

set.seed(420)  

n <- nrow(insurance_data)
train_size <- round(0.6 * n)
valid_size <- round(0.2 * n)
insurance_data_shuffled <- insurance_data[sample(nrow(insurance_data)), ]
train_data <- insurance_data_shuffled[1:train_size, ]
summary(data)
summary(train_data)
valid_data <- insurance_data_shuffled[(train_size + 1):(train_size + valid_size), ]
test_data <- insurance_data_shuffled[(train_size + valid_size + 1):n, ]

# Training tmb_nb on train_data
tmb_nb_train <- glmmTMB(Face ~s(Income_log)+s(TotIncome_log),family =  tweedie() , data= train_data)

xi=family_params(tmb_nb_train)

# Training gamm_nb on train_data
gamm_nb_train <- gamm4(Face~(Income_log)+s(TotIncome_log),family=Tweedie(p=xi), 
                       data = train_data, REML=FALSE) 
summary(tmb_nb_train)
summary(gamm_nb_train$mer)

AIC(tmb_nb_train)
AIC(gamm_nb_train$mer)



set.seed(123)
n_folds <- 5
folds <- cut(seq(1, nrow(valid_data)), breaks = n_folds, labels = FALSE)

# Initialize lists for RMSE
rmse_gamm <- vector("list", n_folds)
rmse_glmmTMB <- vector("list", n_folds)

# Loop through each fold
for(i in 1:n_folds) {
  # Define training and test sets for this fold
  test_indices <- which(folds == i)
  train_indices <- setdiff(seq_len(nrow(valid_data)), test_indices)
  train_fold <- valid_data[train_indices, ]
  test_fold <- valid_data[test_indices, ]
  
  # Train GAMM model
  model_gamm <- gamm4(Face ~s(Income_log)+s(TotIncome_log), family=Tweedie(p=xi),
                      data = train_fold)
  pred_gamm <- predict(model_gamm$gam, test_fold, type = "response")
  rmse_gamm[[i]] <- sqrt(mean((pred_gamm - test_fold$Face)^2))
  
  
  # Train GLMMTMB model
  model_glmmTMB <- glmmTMB(Face ~s(Income_log)+s(TotIncome_log), family =  tweedie(),
                           data = train_fold)
  pred_glmmTMB <- predict(model_glmmTMB, test_fold, type = "response")
  rmse_glmmTMB[[i]] <- sqrt(mean((pred_glmmTMB - test_fold$Face)^2))
}

# Calculate the average RMSE for each model
mean_rmse_gamm <- mean(unlist(rmse_gamm))
mean_rmse_glmmTMB <- mean(unlist(rmse_glmmTMB))

# Print the results
cat("Average RMSE for GAMM: ", mean_rmse_gamm, "\n")
cat("Average RMSE for GLMMTMB: ", mean_rmse_glmmTMB, "\n")


full_train_data <- rbind(train_data, valid_data)

start_time_gamm_train <- Sys.time()
full_model_gamm <- gamm4(Face~s(Income_log)+s(TotIncome_log),family=Tweedie(p=xi), 
                         data = full_train_data)

end_time_gamm_train <- Sys.time()
train_duration_gamm <- end_time_gamm_train - start_time_gamm_train
start_time_gamm_pred <- Sys.time()
predictions_gamm <- predict(full_model_gamm$gam, test_data, type = "response")
end_time_gamm_pred <- Sys.time()
pred_duration_gamm <- end_time_gamm_pred - start_time_gamm_pred


fixed_effects <- fixef(tmb_nb_train)
conditional_beta <- fixed_effects[["cond"]]
#random_effects <- ranef(tmb_nb_train)

start_time_glmmTMB_train <- Sys.time()

full_model_glmmTMB <- glmmTMB(Face~s(Income_log)+s(TotIncome_log),tweedie(), 
                              data = full_train_data,REML=TRUE)

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

TMBRMSE <- sqrt(mean((predictions_glmmTMB - test_data$Face)^2))
GAMMRMSE<-sqrt(mean((predictions_gamm - test_data$Face)^2))

print(paste("RMSE glmmTMB:", TMBRMSE))
print(paste("RMSE gamm4:", GAMMRMSE))




# Correctly combine predictions for the plot
plot_data <- data.frame(
  Actual = rep(test_data$Face, 2),
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


actual_density <- density(full_train_data$Face)

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


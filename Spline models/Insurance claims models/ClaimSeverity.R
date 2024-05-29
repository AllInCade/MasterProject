source("requirements.R")

data("motorins")
view(motorins)
?motorins
head(motorins)

hist(motorins$Claims, breaks=2000)

motorins$HasClaim <- ifelse(motorins$Claims > 0, 1, 0)
n <- nrow(motorins)
train_size <- round(0.6 * n)
valid_size <- round(0.2 * n)
data_shuffled <- motorins[sample(nrow(motorins)), ]
train_data <- data_shuffled[1:train_size, ]
valid_data <- data_shuffled[(train_size + 1):(train_size + valid_size), ]
test_data <- data_shuffled[(train_size + valid_size + 1):n, ]


classifier <- glmmTMB(HasClaim ~ s(Insured)+(1|Zone),
                                        family = binomial(link = "logit"), data = train_data)

summary(classifier)
valid_data$predicted_class <- predict(classifier, newdata = valid_data, type = "response", allow.new.levels=TRUE)
valid_data$predicted_class <- ifelse(valid_data$predicted_class >= 0.5, 1, 0)

confusionMatrix(as.factor(valid_data$predicted_class), as.factor(valid_data$HasClaim))





motorins<-motorins[motorins$Claims>0,]

mean(motorins$Payment)

# Finding most important variables

set.seed(101) 
training_indices <- sample(1:nrow(motorins), 0.8 * nrow(motorins))
train_data <- motorins[training_indices, ]
test_data <- motorins[-training_indices, ]

model <- randomForest(Payment~Insured+Make+Bonus+Kilometres, data = train_data, na.action = na.omit, importance = TRUE)

predictions <- predict(model, test_data)
actuals <- test_data$Payment



RMSE <- sqrt(mean((predictions - actuals)^2))
print(paste("RMSE:", RMSE))



importance_scores <- randomForest::importance(model, type = 2) 
feature_names <- row.names(importance_scores)
importance_data <- data.frame(Feature = feature_names, Importance = importance_scores[, 1])
sorted_importance_data <- importance_data[order(-importance_data$Importance), ] # Sort in decreasing order
barplot(sorted_importance_data$Importance, names.arg = sorted_importance_data$Feature, 
        main = "Feature Importance in Predicting mdeath", las = 2, cex.names = 0.7)

print(sorted_importance_data)





glmmTMB_model<-glmmTMB(Payment~s(Insured)+Bonus+(1|Zone)+offset(log(Claims)),data=motorins,family = Gamma(link="log"))



gamm4_model<- gamm4(Payment~s(Insured)+Bonus,random=~(1|Zone)+offset(log(Claims)),data=motorins,family=Gamma(link="log"))





summary(glmmTMB_model)

summary(gamm4_model$mer)

#As per now the models of gamm4 and glmmtmb are identical



# Data split 60/20/20

set.seed(420)  

n <- nrow(motorins)
train_size <- round(0.6 * n)
valid_size <- round(0.2 * n)
motorins_shuffled <- motorins[sample(nrow(motorins)), ]
train_data <- motorins_shuffled[1:train_size, ]
valid_data <- motorins_shuffled[(train_size + 1):(train_size + valid_size), ]
test_data <- motorins_shuffled[(train_size + valid_size + 1):n, ]

# Training tmb_nb on train_data
tmb_nb_train <- glmmTMB(Payment~s(Insured)+Bonus+ (1|Zone)+offset(log(Claims)),disp=~Kilometres,family = Gamma(link = "log"),
                        data = train_data)


# Training gamm_nb on train_data
gamm_nb_train <- gamm4(Payment~s(Insured)+Bonus+offset(log(Claims)),random=~(1|Zone),family = Gamma(link = "log"), 
                       data = train_data) 

summary(tmb_nb_train)
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
  model_gamm <- gamm4(Payment~s(Insured)+Bonus+offset(log(Claims)),family = Gamma(link = "log"), 
                      data = train_fold)
  pred_gamm <- predict(model_gamm$gam, test_fold, type = "response")
  rmse_gamm[[i]] <- sqrt(mean((pred_gamm - test_fold$Payment)^2))
  
  # Train GLMMTMB model
  model_glmmTMB <- glmmTMB(Payment~s(Insured)+Bonus+offset(log(Claims)),disp=~Kilometres,family = Gamma(link = "log"),
                           data = train_fold)
  pred_glmmTMB <- predict(model_glmmTMB, test_fold, type = "response")
  rmse_glmmTMB[[i]] <- sqrt(mean((pred_glmmTMB - test_fold$Payment)^2))
}

# Calculate the average RMSE for each model
mean_rmse_gamm <- mean(unlist(rmse_gamm))
mean_rmse_glmmTMB <- mean(unlist(rmse_glmmTMB))

# Print the results
cat("Average RMSE for GAMM: ", mean_rmse_gamm, "\n")
cat("Average RMSE for GLMMTMB: ", mean_rmse_glmmTMB, "\n")


full_train_data <- rbind(train_data, valid_data)

start_time_gamm_train <- Sys.time()


full_model_gamm <- gamm4(Payment~s(Insured)+Bonus+offset(log(Claims)),family = Gamma(link="log"), data = full_train_data)

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


hist(motorins$Payment, breaks = 2000)

full_model_glmmTMB <- glmmTMB(Payment~s(Insured)+Bonus+offset(log(Claims)),disp=~Kilometres,family = Gamma(link="log"), 
                              data = full_train_data)


summary(full_model_glmmTMB)

plot(resid(full_model_glmmTMB))

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

AIC(full_model_glmmTMB)
AIC(full_model_gamm$mer)
actuals <- test_data$Payment

RMSETMB <- sqrt(mean((predictions_glmmTMB - actuals)^2))

RMSEGAMM <- sqrt(mean((predictions_gamm - actuals)^2))

print(paste("RMSE for glmmTMB:", RMSETMB))

print(paste("RMSE for gamm4:", RMSEGAMM))

# Correctly combine predictions for the plot
plot_data <- data.frame(
  Actual = rep(test_data$Payment, 2),
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


actual_density <- density(full_train_data$Payment)

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
  facet_wrap(~Model)+
  coord_cartesian(xlim = c(0,1000000)) 


source("requirements.R")

swedish_mortality<-read_sas('mortality.sas7bdat')

head(swedish_mortality)

summary(swedish_mortality)


# Assuming swedish_mortality is your dataframe and Male_death is the variable
# Round to the nearest integer, then convert to integer
swedish_mortality$Male_death <- as.integer(round(swedish_mortality$Male_death))

swedish_mortality$Male_Exp <- as.integer(round(swedish_mortality$Male_Exp))
# Assuming 'data' is your dataset and it includes a column 'death_rate_q_male'
# that represents the mortality rate for males, you can create a subset of
# the data where the death rate is strictly between 0 and 1.

swedish_mortality1 <- swedish_mortality[swedish_mortality$q_male > 0 & swedish_mortality$q_male < 1,]


na_counts <- sapply(swedish_mortality1, function(x) sum(is.na(x)))

# Print the counts of NA values
print(na_counts)

swedish_mortality1<-na.omit(swedish_mortality1)

swedish_mortality1$obs <- factor(seq(nrow(swedish_mortality1)))



set.seed(101) 
training_indices <- sample(1:nrow(swedish_mortality1), 0.8 * nrow(swedish_mortality1))
train_data <- swedish_mortality1[training_indices, ]
test_data <- swedish_mortality1[-training_indices, ]

model <- randomForest(q_male ~  Year+Age+Female_death+ Male_Exp+Female_Exp+ Male_death+q_female, data = train_data, na.action = na.omit, importance = TRUE)

predictions <- predict(model, test_data)
actuals <- test_data$q_male



RMSE <- sqrt(mean((predictions - actuals)^2))
print(paste("RMSE:", RMSE))



importance_scores <- randomForest::importance(model, type = 2) 
feature_names <- row.names(importance_scores)
importance_data <- data.frame(Feature = feature_names, Importance = importance_scores[, 1])
sorted_importance_data <- importance_data[order(-importance_data$Importance), ] # Sort in decreasing order
barplot(sorted_importance_data$Importance, names.arg = sorted_importance_data$Feature, 
        main = "Feature Importance in Predicting mdeath", las = 2, cex.names = 0.7)

print(sorted_importance_data)



# Assuming df is your dataframe and male_death_count is your variable
num_bins <- 5 # For example, creating 5 bins

# Create quantile bins and convert to factor
swedish_mortality1$death_count_group <- cut(swedish_mortality1$Male_death, 
                            breaks = quantile(swedish_mortality1$Male_death, probs = seq(0, 1, length.out = num_bins + 1), na.rm = TRUE),
                            include.lowest = TRUE,
                            labels = FALSE)

# Convert numeric bins to factor for use as a random effect
swedish_mortality1$death_count_group <- as.factor(swedish_mortality1$death_count_group)


swedish_mortality1$Age_scaled <- scale(swedish_mortality1$Age)




swedish_mortality1$q_male_logit <- log(swedish_mortality1$q_male / (1 - swedish_mortality1$q_male))
hist(swedish_mortality1$q_male_logit)


gamm_model <- gamm4(q_male_logit~ s(Age) +s(Year),random=~(1|death_count_group),
                   family=gaussian(),
                   data = swedish_mortality1) 


glmmTMB_model <- glmmTMB(q_male ~ s(Year) + s(Age)+ (1|death_count_group),
                         family=beta_family(link = "logit"), 
                         data = swedish_mortality1)


# Check the summary
summary(gamm_model$mer)

summary(glmmTMB_model)






# Data split 60/20/20

set.seed(420)  

n <- nrow(swedish_mortality1)
train_size <- round(0.6 * n)
valid_size <- round(0.2 * n)
mortality_data_shuffled <- swedish_mortality1[sample(nrow(swedish_mortality1)), ]
train_data <- mortality_data_shuffled[1:train_size, ]
valid_data <- mortality_data_shuffled[(train_size + 1):(train_size + valid_size), ]
test_data <- mortality_data_shuffled[(train_size + valid_size + 1):n, ]

# Training tmb_nb on train_data
tmb_nb_train <- glmmTMB(q_male ~ s(Year)+s(Age)+ (1|death_count_group),
                        family=beta_family(link="logit"), 
                        data = train_data)


# Training gamm_nb on train_data
gamm_nb_train <- gamm4(q_male_logit~ s(Age) +s(Year),random=~(1|death_count_group),
                       family=gaussian(),
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
  model_gamm <- gamm4(q_male_logit~ s(Age) +s(Year),random=~(1|death_count_group),
                      family=gaussian(),
                      data = train_fold, REML=TRUE)
  # Make predictions on the logit scale
  logit_predictions <- predict(model_gamm$gam, test_fold, type = "response")
  # Back-transform predictions to the original scale
  pred_gamm  <- 1 / (1 + exp(-logit_predictions))
  rmse_gamm[[i]] <- sqrt(mean((pred_gamm - test_fold$q_male)^2))
  
  # Train GLMMTMB model
  model_glmmTMB <- glmmTMB(q_male~ s(Age)+s(Year)+ (1|death_count_group), family=beta_family(link = "logit"), 
                           data = train_fold, REML=TRUE)
  pred_glmmTMB <- predict(model_glmmTMB, test_fold, type = "response")
  rmse_glmmTMB[[i]] <- sqrt(mean((pred_glmmTMB - test_fold$q_male)^2))
}




# Calculate the average RMSE for each model
mean_rmse_gamm <- mean(unlist(rmse_gamm))
mean_rmse_glmmTMB <- mean(unlist(rmse_glmmTMB))

# Print the results
cat("Average RMSE for GAMM: ", mean_rmse_gamm, "\n")
cat("Average RMSE for GLMMTMB: ", mean_rmse_glmmTMB, "\n")


full_train_data <- rbind(train_data, valid_data)

start_time_gamm_train <- Sys.time()
full_model_gamm <- gamm4(q_male_logit~ s(Age) +s(Year),random=~(1|death_count_group),
                         family=gaussian(),
                         data = full_train_data)




end_time_gamm_train <- Sys.time()
train_duration_gamm <- end_time_gamm_train - start_time_gamm_train
start_time_gamm_pred <- Sys.time()
# Make predictions on the logit scale
logit_predictions <- predict(full_model_gamm$gam, test_data, type = "response")
# Back-transform predictions to the original scale
predictions_gamm  <- 1 / (1 + exp(-logit_predictions))
end_time_gamm_pred <- Sys.time()
pred_duration_gamm <- end_time_gamm_pred - start_time_gamm_pred


fixed_effects <- fixef(tmb_nb_train)
conditional_beta <- fixed_effects[["cond"]]
#random_effects <- ranef(tmb_nb_train)

start_time_glmmTMB_train <- Sys.time()


full_model_glmmTMB <- glmmTMB(q_male~ s(Age) +s(Year)+ (1|death_count_group),
                              family=beta_family(link = "logit"), 
                              data = full_train_data)

end_time_glmmTMB_train <- Sys.time()
train_duration_glmmTMB <- end_time_glmmTMB_train - start_time_glmmTMB_train

start_time_glmmTMB_pred <- Sys.time()
predictions_glmmTMB <- predict(full_model_glmmTMB, test_data, type = "response", allow.new.levels=TRUE)
end_time_glmmTMB_pred <- Sys.time()
pred_duration_glmmTMB <- end_time_glmmTMB_pred - start_time_glmmTMB_pred
cat("Training duration for GAM: ", train_duration_gamm, "\n")
cat("Training duration for GLMMTMB: ", train_duration_glmmTMB, "\n")

cat("Prediction duration for GAM: ", pred_duration_gamm, "\n")
cat("Prediction duration for GLMMTMB: ", pred_duration_glmmTMB, "\n")



AIC(full_model_glmmTMB)
AIC(full_model_gamm$mer)
actuals <- test_data$q_male



RMSETMB <- sqrt(mean((predictions_glmmTMB - actuals)^2))

RMSEGAMM <- sqrt(mean((predictions_gamm - actuals)^2))

print(paste("RMSE for glmmTMB:", RMSETMB))
print(paste("RMSE for gam:", RMSEGAMM))


# Correctly combine predictions for the plot
plot_data <- data.frame(
  Actual = rep(test_data$q_male, 2),
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


actual_density <- density(full_train_data$q_male)

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




library(haven)
library(glmmTMB)
library(mgcv)
library(ggplot2)
library(dplyr)
library(DHARMa)
library(corrplot)
library(tidyverse)
library(readr)
library(tidyr)
library(caret)
library(AICcmodavg)
library(randomForest)
library(importance)
library(gamm4)

swedish_mortality<-read_sas('mortality.sas7bdat')

head(swedish_mortality)

na_counts <- sapply(swedish_mortality, function(x) sum(is.na(x)))

# Print the counts of NA values
print(na_counts)

# Key Variables including Education
key_vars <- c("Year"  , "Age" ,"Female_Exp" ,"Male_Exp", "q_female",   "q_male", "Female_death" ,"Male_death", "L_female_exp", "L_male_exp")

# Correlation Matrix including Education
cov_matrix_extended <- cor(swedish_mortality[key_vars])

# Plotting the covariance matrix
cov_plot <- corrplot(cov_matrix_extended, method = "number", type = "upper", tl.col = "black")

# Round Male_death to the nearest integer
swedish_mortality$Male_death <- round(swedish_mortality$Male_death)


# Assuming swedish_mortality is your dataframe and Male_death is the variable
# Round to the nearest integer, then convert to integer
swedish_mortality$Male_death <- as.integer(round(swedish_mortality$Male_death))

swedish_mortality<-na.omit(swedish_mortality)
# 3. Visualization
# Age Distribution
g1<-ggplot(data = swedish_mortality, aes(x = Age)) + 
  geom_histogram(binwidth = 1, fill = "blue", color = "black") +
  labs(title = "Age Distribution", x = "Age", y = "Count")

# Income Distribution
g2<-ggplot(data = swedish_mortality, aes(x = Year)) + 
  geom_histogram(binwidth = 1, fill = "green", color = "black") +
  labs(title = "Year Distribution", x = "Year", y = "Count")

# Face Value Distribution
g3<-ggplot(data = swedish_mortality, aes(x = Female_death)) + 
  geom_histogram(binwidth = 1, fill = "red", color = "black") +
  labs(title = "Female_Death Distribution", x = "Female_Death", y = "Count")

# Face Value Distribution
g4<-ggplot(data = swedish_mortality, aes(x = Male_death)) + 
  geom_histogram(binwidth = 1, fill = "red", color = "black") +
  labs(title = "male_Death Distribution", x = "Male_death", y = "Count")


gridExtra::grid.arrange(g1, g2, g3,g4, ncol = 4)




# Finding most important variables

set.seed(101) 
training_indices <- sample(1:nrow(swedish_mortality), 0.8 * nrow(swedish_mortality))
train_data <- swedish_mortality[training_indices, ]
test_data <- swedish_mortality[-training_indices, ]

model <- randomForest(Male_death ~  Year+Age+Female_death+ Male_Exp+Female_Exp+ q_male+q_female, data = train_data, na.action = na.omit, importance = TRUE)

predictions <- predict(model, test_data)
actuals <- test_data$Male_death



RMSE <- sqrt(mean((predictions - actuals)^2))
print(paste("RMSE:", RMSE))



importance_scores <- randomForest::importance(model, type = 2) 
feature_names <- row.names(importance_scores)
importance_data <- data.frame(Feature = feature_names, Importance = importance_scores[, 1])
sorted_importance_data <- importance_data[order(-importance_data$Importance), ] # Sort in decreasing order
barplot(sorted_importance_data$Importance, names.arg = sorted_importance_data$Feature, 
        main = "Feature Importance in Predicting mdeath", las = 2, cex.names = 0.7)

print(sorted_importance_data)



poisson_model <- glm(Male_death ~ Age + Year+ log(Male_Exp) ,family = poisson(),
                     data = swedish_mortality)
# Summary of the model
summary(poisson_model)

library(AER)

dispersiontest(poisson_model)

#Testing if negative binomial model will work better

nbin_model<-glm(Male_death ~ Age + Year+log(Male_Exp),family = negbin(theta = 0.773),
                data = swedish_mortality)

summary(nbin_model)


glmmTMB_model<- glmmTMB(Male_death~ s(Age) +s(Year), offset=log(Male_Exp),family = nbinom2(), 
                        data = swedish_mortality)

gamm4_model<- gamm4(Male_death~ s(Age) +s(Year)+offset(log(Male_Exp)),family = negbin(theta=10.2),
                    data = swedish_mortality)


summary(glmmTMB_model)

summary(gamm4_model$mer)

#As per now the models of gamm4 and glmmtmb are identical





# Data split 60/20/20

set.seed(420)  

n <- nrow(swedish_mortality)
train_size <- round(0.6 * n)
valid_size <- round(0.2 * n)
mortality_data_shuffled <- swedish_mortality[sample(nrow(swedish_mortality)), ]
train_data <- mortality_data_shuffled[1:train_size, ]
valid_data <- mortality_data_shuffled[(train_size + 1):(train_size + valid_size), ]
test_data <- mortality_data_shuffled[(train_size + valid_size + 1):n, ]

# Training tmb_nb on train_data
tmb_nb_train <- glmmTMB(Male_death~ s(Age) +s(Year),offset=log(Male_Exp),disp=~Age,family = nbinom2(), 
                        data = train_data)


# Training gamm_nb on train_data
gamm_nb_train <- gamm4(Male_death~ s(Age) +s(Year)+offset(log(Male_Exp)),family = negbin(theta=10.2), 
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
  model_gamm <- gamm4(Male_death~ s(Age) +s(Year)+ offset(log(Male_Exp)),family = negbin(theta=10.2), 
                      data = train_fold)
  pred_gamm <- predict(model_gamm$gam, test_data, type = "response")
  rmse_gamm[[i]] <- sqrt(mean((pred_gamm - test_data$Male_death)^2))
  
  # Train GLMMTMB model
  model_glmmTMB <- glmmTMB(Male_death~ s(Age) +s(Year),offset=log(Male_Exp),disp=~Age,family = nbinom2(),
                           data = train_fold)
  pred_glmmTMB <- predict(model_glmmTMB, test_data, type = "response")
  rmse_glmmTMB[[i]] <- sqrt(mean((pred_glmmTMB - test_data$Male_death)^2))
}

# Calculate the average RMSE for each model
mean_rmse_gamm <- mean(unlist(rmse_gamm))
mean_rmse_glmmTMB <- mean(unlist(rmse_glmmTMB))

# Print the results
cat("Average RMSE for GAMM: ", mean_rmse_gamm, "\n")
cat("Average RMSE for GLMMTMB: ", mean_rmse_glmmTMB, "\n")


full_train_data <- rbind(train_data, valid_data)

start_time_gamm_train <- Sys.time()
full_model_gamm <- gamm4(Male_death~s(Age) +s(Year)+offset(log(Male_Exp)),family = negbin(theta=10.2), 
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

full_model_glmmTMB <- glmmTMB(Male_death~ s(Age) +s(Year),offset = log(Male_Exp),disp=~Age,family = nbinom2(), 
                              data = full_train_data)

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
actuals <- test_data$Male_death



RMSETMB <- sqrt(mean((predictions_glmmTMB - actuals)^2))

RMSEGAMM <- sqrt(mean((predictions_gamm - actuals)^2))

print(paste("RMSE for glmmTMB:", RMSETMB))
print(paste("RMSE for gamm4:", RMSEGAMM))


# Correctly combine predictions for the plot
plot_data <- data.frame(
  Actual = rep(test_data$Male_death, 2),
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


actual_density <- density(full_train_data$Male_death)

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



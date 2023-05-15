rm(list = ls())
set.seed(123)
library(plotly)
library("pROC")
file_path=file.path("C:/Users/ldegr/Desktop/Thesis/ASML case PBAS 2023/ASML_case_PBAS 2023.csv")
data=read.csv(file_path)
data=data.frame(data)
library(plotly)
library(corrplot)
col_names <- colnames(data)
for (col_name in col_names) {
  assign(col_name, data[[col_name]])
}
yield <- Yieldproxy_y_perc
data$wafer_id <- NULL
wafer_x_microm <- wafer_coordinate_x_mm*1000
wafer_y_microm <- wafer_coordinate_y_mm*1000
fieldcenter_x_microm <- field_center_x
fieldcenter_y_microm <- field_center_y
wafer_x_abs <- wafer_x_microm
field_x <- field_center_x
field_y <- field_center_y
ovl_nm<- abs(OVL_y_nm)
disttocenter_wafer <- sqrt(wafer_x_microm^2 + wafer_y_microm^2)
disttocenter_field <- sqrt(intra_field_x^2 + intra_field_y^2)
fieldcenterdistance <- sqrt(field_center_x^2 + field_center_y^2)
fieldcenterdistancey <- sqrt(intra_field_y^2 + intra_field_y^2)
hist(yield, breaks = 100) 
hist(yield, breaks = 10) 
hist(field_x, breaks = 100)
hist(field_y)
hist(wafer_x_microm)
hist(wafer_y_microm, breaks = 20)
# count where yield < 80
sum(yield < 80, na.rm = TRUE)
sum(yield < 50, na.rm = TRUE)

# plot
fig <- plot_ly(data, x = wafer_x_microm, y = wafer_y_microm, z = yield, type = "scatter3d", mode = "markers")
fig <- fig %>% layout(title = "Yield vs. Wafer Coordinates", scene = list(xaxis = list(title = "Wafer X"), yaxis = list(title = "Wafer Y"), zaxis = list(title = "Yield")))
fig

# correlation plot
#install.packages("corrplot")


# regression
fit1 <- lm(yield ~ ovl_nm+ intra_field_x + L1_CD_y_nm + L2_CD_y_nm)
summary(fit1)

hist(data$field_center_x, breaks =100)
library(dplyr)
n_distinct(data$field_center_x,data$field_center_y)

fit2 <- lm(yield ~ ovl_nm + L1_CD_y_nm + L2_CD_y_nm +disttocenter_field+disttocenter_wafer)
summary(fit2)
#create hardcopy drop yield
name= names(data) %in% c("Design_y_nm","intra_field_x","intra_field_y","wafer_coordinate_x_mm","wafer_coordinate_y_mm","field_center_x","field_center_y") 
data_copy=data[!name]
data_copy['field_center_x_microm']=fieldcenter_x_microm
data_copy['field_center_y_microm']=fieldcenter_y_microm
data_copy['absovly']=ovl_nm
data.cor = cor(data_copy)
corrplot(data.cor)
#test en train set maken
train.size=dim(data_copy)[1]/2
train=sample(1:dim(data_copy)[1],train.size)
test=-train
data_copy.train=data_copy[train,]
data_copy.test=data_copy[test,]
train.mat=model.matrix(Yieldproxy_y_perc~.,data=data_copy.train)
test.mat=model.matrix(Yieldproxy_y_perc~.,data=data_copy.test)

library(glmnet)
name= names(data) %in% c("Yieldproxy_y_perc","Design_y_nm","intra_field_x","intra_field_y","wafer_coordinate_x_mm","wafer_coordinate_y_mm","field_center_x","field_center_y","OVL_y_nm") 
data_copy2=data[!name]
data_copy2['field_center_x_microm']=fieldcenter_x_microm
data_copy2['field_center_y_microm']=fieldcenter_y_microm
data_copy2['absovly']=ovl_nm
data_copy3=as.matrix(data_copy2)
fit3=cv.glmnet(data_copy3, yield, alpha=.5)
best_lambda <- fit3$lambda.min
best_lambda
plot(fit3) 

best_model <- glmnet(data_copy2, yield, alpha=.5, lambda = best_lambda)
coef(best_model)

data['BinVar']=ifelse(data['Yieldproxy_y_perc']<98,0,1) #HIER PERCENTAGE AANPASSEN
name= names(data) %in% c("Yieldproxy_y_perc","Design_y_nm") 
data_copy4=data[!name] #verwijder yield en nutteloos 
data_copy4[['BinVar']]=as.matrix(data_copy4['BinVar'])

data_copy4$BinVar <- sapply(data_copy4$BinVar, unlist) 
data_copy4$BinVar <- as.double(data_copy4$BinVar)
data_copy_4subset <- data_copy4[, !(colnames(data_copy4) %in% "BinVar")]

data_copy_4subset <- as.matrix(data_copy_4subset, 1470,9)  #HIER KPI AAN TOEVOEGEN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

data_copy_4subset <- apply(data_copy_4subset, 2, as.numeric)
hist(data_copy4$BinVar, 2)  


# Extract the response variable as a vector
# Convert data frame to matrix


# Fit the Lasso logistic regression model
fit_cv <- cv.glmnet(data_copy_4subset , data_copy4$BinVar, family = "binomial", alpha =0 )
best_lambda <- fit_cv$lambda.min
bestfit=glmnet(data_copy_4subset , data_copy4$BinVar, family = "binomial", alpha = 0, lambda=best_lambda)
# Fit the Lasso logistic regression model
coef(bestfit)




##Beste versie





data_copy_4subset2=data.frame(data_copy_4subset)
#data_copy_4subset2$L1_over_L2 <- data_copy_4subset2$L1_CD_y_nm / data_copy_4subset2$L2_CD_y_nm
#data_copy_4subset2$KPI=abs(data_copy_4subset2$L1_CD_y_nm - data_copy_4subset2$L2_CD_y_nm)*abs(data_copy_4subset2$OVL_y_nm)
data_copy_4subset2$absOVL_y_nm=abs(data_copy_4subset2$OVL_y_nm)
# Exclude the two original columns
data_copy_4subset2$MaxoverMinL1_over_L2 <- pmax(data_copy_4subset2$L1_CD_y_nm,data_copy_4subset2$L2_CD_y_nm)/ pmin(data_copy_4subset2$L1_CD_y_nm,data_copy_4subset2$L2_CD_y_nm)
data_copy_4subset2$maxminabs=data_copy_4subset2$MaxoverMinL1_over_L2*data_copy_4subset2$absOVL_y_nm
data_copy_4subset2$disttocenter_fieldabsoverlay=data_copy_4subset2$absOVL_y_nm*disttocenter_field
data_copy_4subset2$disttocenter_fielddmaxmin=data_copy_4subset2$MaxoverMinL1_over_L2 *disttocenter_field
data_copy_4subset2$fieldcenterdistanceabsoverlay=data_copy_4subset2$absOVL_y_nm*fieldcenterdistance
data_copy_4subset2$fieldcenterdistancecdmaxmin=data_copy_4subset2$MaxoverMinL1_over_L2 **fieldcenterdistance
data_copy_4subset2$fieldcenter=fieldcenterdistance
data_copy_4subset2$disttocenter=disttocenter_field
data_copy_5subset = data_copy_4subset2[ , !(names(data_copy_4subset2) %in% c("L1_CD_y_nm", "L2_CD_y_nm","OVL_y_nm","wafer_coordinate_x_mm","wafer_coordinate_y_mm", "field_center_x", "field_center_y", "intra_field_x", "intra_field_y"))]


# Convert to matrix
data_copy5subset_matrix = as.matrix(data_copy_5subset)



library(caret)
library(pROC)

preprocessed_data <- preProcess(data_copy5subset_matrix, method = c("center", "scale"))
scaled_data <- predict(preprocessed_data, data_copy5subset_matrix)

data_copy4$BinVarTwoFactor <- factor(data_copy4$BinVar, levels = c(0, 1), labels = c("No", "Yes"))
# Convert outcome variable to a two-level factor
data_copy4$BinVarTwoFactor <- factor(ifelse(data_copy4$BinVar == 0, "No", "Yes"))

# Make predictions with the Bagging model


# Set up train control with bootstrapping
library(caret)
library(pROC)

# ... [previous code omitted for brevity] ...

f1Summary <- function(data, lev = NULL, model = NULL) {
  cm <- confusionMatrix(data$pred, data$obs, mode = "prec_recall")
  f1_val <- 2 * (cm$byClass["Precision"] * cm$byClass["Recall"]) / (cm$byClass["Precision"] + cm$byClass["Recall"])
  out <- c(F1 = f1_val)
  names(out) <- "F1"
  out
}

boot_control <- trainControl(
  method = "boot",
  number = 100,
  savePredictions = "all",
  classProbs = TRUE,
  summaryFunction = f1Summary
)

logit_ridge_model <- train(
  x = scaled_data,
  y = data_copy4$BinVarTwoFactor,
  method = "glmnet",
  metric = "F1",
  trControl = boot_control,
  family = "binomial",
  tuneGrid = expand.grid(alpha = .5, lambda = seq(0.01, 0.1, length.out = 10))
)

# Make predictions with the logistic regression model
predictions_ridge <- predict(logit_ridge_model, scaled_data, type = "prob")

# Calculate AUC-ROC for Ridge
roc_obj_ridge <- roc(data_copy4$BinVarTwoFactor, predictions_ridge[, "Yes"], levels = rev(levels(data_copy4$BinVarTwoFactor)))
auc_roc_ridge <- as.numeric(auc(roc_obj_ridge))



# Get resampled predictions for each cutoff value
cutoff_values <- seq(0.1, 0.99, by = 0.01)

f1_scores_ridge <- vector("numeric", length(cutoff_values))

for (i in 1:length(cutoff_values)) {
  cutoff <- cutoff_values[i]
  resample_f1_scores_ridge <- vector("numeric", 100)
  
  for (j in 1:100) {
    resample_name <- sprintf("Resample%02d", j)
    
    resample_predictions_ridge <- logit_ridge_model$pred[logit_ridge_model$pred$Resample == resample_name, ]
    
    
    y_true_ridge <- factor(resample_predictions_ridge$obs, levels = c("No", "Yes"))
    
    # Bagging model F1 score
    
    
    # Ridge model F1 score
    y_pred_ridge <- factor(ifelse(resample_predictions_ridge$Yes >= cutoff, "Yes", "No"), levels = c("No", "Yes"))
    cm_ridge <- table(Predicted = y_pred_ridge, Actual = y_true_ridge)
    precision_ridge <- cm_ridge[2, 2] / (sum(cm_ridge[, 2]))
    recall_ridge <- cm_ridge[2, 2] / (sum(cm_ridge[2, ]))
    resample_f1_scores_ridge[j] <- 2 * (precision_ridge * recall_ridge) / (precision_ridge + recall_ridge)
  }
  
  
  f1_scores_ridge[i] <- mean(resample_f1_scores_ridge, na.rm = TRUE)
  
  cat(cutoff, "Ridge Average F1 score:", f1_scores_ridge[i], "\n")
}

# ... continue with the rest of the code ...

# Find the best cutoff and F1 score for Bagging



# Find the best cutoff and F1 score for Ridge
best_index_ridge <- which.max(f1_scores_ridge)
best_cutoff_ridge <- cutoff_values[best_index_ridge]
best_avg_f1_score_ridge <- f1_scores_ridge[best_index_ridge]

cat("Best cutoff for Elastic Net:", best_cutoff_ridge, "with the highest average F1 score:", best_avg_f1_score_ridge, "\n")


est_lambda_ridge <- logit_ridge_model$bestTune$lambda

# Train the final Ridge model with the optimal lambda value on the entire dataset
final_ridge_model <- train(
  x = scaled_data,
  y = data_copy4$BinVarTwoFactor,
  method = "glmnet",
  metric = "ROC",
  trControl = trainControl(method = "none", classProbs = TRUE),
  family = "binomial",
  tuneGrid = expand.grid(alpha = .5, lambda = est_lambda_ridge)
)

# Get the coefficients from the final Ridge model
final_ridge_coefs <- coef(final_ridge_model$finalModel, s = est_lambda_ridge)
final_ridge_coefs


data_5_subset2_cor <- cor(data_copy5subset_matrix)

corrplot(data_5_subset2_cor)
 

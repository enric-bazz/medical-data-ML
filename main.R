# Clear command
cat("\014")
# Clear environment
rm(list=ls())
#Setting working directory
dir <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(dir))

## Loading data
load('data_final_project.Rdata')
cat('Number of observations =', nrow(data))
cat('Number of variables =', ncol(data))
summary(data) #shows NAs in each column

table(data$outcome) #checks for class unbalance
cat('Fraction of observations with outcome 0 =', table(data$outcome)[1]/nrow(data))
cat('Fraction of observations with outcome 1 =', table(data$outcome)[2]/nrow(data))

### Training/Test set split
ind_0 <- which(data$outcome==0) #indexes of 0 class
ind_1 <- which(data$outcome==1) #indexes of 1 class
s <- 0.8 #splitting fraction

set.seed(123) #setting RNG seed for reproducibility
ind_0_training <- sample(ind_0, round(s*table(data$outcome)[1]), replace=F) #sampling from ind_0
ind_1_training <- sample(ind_1, round(s*table(data$outcome)[2]), replace=F) #sampling from ind_1
ind_training <- c(ind_0_training,ind_1_training) #concatenate training set indexes
training_set <- data[ind_training,]
cat('Fraction of observations with outcome 0 in training =', table(training_set$outcome)[1]/nrow(training_set))
cat('Fraction of observations with outcome 1 =', table(training_set$outcome)[2]/nrow(training_set))
test_set <- data[-ind_training,] 
cat('Fraction of observations with outcome 0 =', table(test_set$outcome)[1]/nrow(test_set))
cat('Fraction of observations with outcome 1 =', table(test_set$outcome)[2]/nrow(test_set))

### Data Pre-processing 
numerical_var <- sapply(training_set, is.numeric) #numerical variables names (numeric and integer)
factor_var <- sapply(training_set, is.factor) #factor variables names
factor_var <- !numerical_var #equivalent

# Extreme data points
for(i in which(numerical_var)){  #conditional indexing of the indices
  boxplot(training_set[,i], main=paste0('Boxplot of ', colnames(training_set)[i]))
}
outliers_var <- c('ferritin','crprot','trig','gluc') #possible outliers
par(mfrow=c(2,2))
for(j in which(colnames(training_set) %in% outliers_var)){
  boxplot(training_set[,j], main=paste0('Boxplot of ', colnames(training_set)[j]))
}
#Remove outliers
outliers_var <- c('crprot','trig','gluc') #definitive outliers
outlier_thresholds <- c(80,300,1000) #manually selected thresholds
ind_outliers <- numeric() #initialize vector
par(mfrow=c(1,3))
for(k in 1:length(outlier_thresholds)){
  ind_var <- which(colnames(training_set) %in% outliers_var)[k] #variable column
  ind_outliers <- c(ind_outliers,which(training_set[,ind_var]>outlier_thresholds[k])) #outlier row (append each time)
  boxplot(training_set[,ind_var],main=paste0('Boxplot of ', colnames(training_set)[ind_var]))
  abline(h=outlier_thresholds[k],col='red')
}

training_set <- training_set[-unique(ind_outliers),]  #eliminate whole observation in training set
#unique() manages repeating of outliers
nrow(training_set) #lower than before
for(i in which(colnames(training_set) %in% outliers_var)){  #check effectiveness of removal
  boxplot(training_set[,i], main=paste0(colnames(training_set)[i]))
}


#Update indices for later use
ind_0_training <- which(training_set$outcome==0)
ind_1_training <- which(training_set$outcome==1)

# Collinearity
par(mfrow=c(1,1))
corr.matrix <- cor(training_set[,numerical_var], use = 'complete.obs') #correlation matrix (using conditional indexing)
#"use = 'complete.obs'" excludes observation with NAs, equivalent to na.omit()

corrplot::corrplot(corr.matrix, method = 'color')
corr_threshold <- 0.8 #correlation coefficient threshold
high_corr <- which(abs(corr.matrix>0.8) & abs(corr.matrix)<1, arr.ind = T) 
#use abs() to include negative correlations

par(mfrow=c(1,2))
graphics::plot(training_set$waist, training_set$bmi, main = 'Scatter plot', xlab = 'waist', ylab = 'bmi') #waist vs bmi
graphics::plot(training_set$ldl, training_set$tot_chol, main = 'Scatter plot', xlab = 'ldl', ylab = 'tot_chol') #ldl vs tot_chol

#Strong but below threshold correlation
par(mfrow=c(1,1))
graphics::plot(training_set$diastolic_bp, training_set$systolic_bp, main = 'Scatter plot', xlab = 'diastolic_bp', ylab = 'systolic_bp')

#Eliminate the correlated variables waist and tot_chol in both sets
training_set$waist <- NULL
test_set$waist <- NULL
training_set$tot_chol <- NULL
test_set$tot_chol <- NULL
# Check new number of features
ncol(training_set)
ncol(test_set)
# Update variables of indices
numerical_var <- sapply(training_set, is.numeric) #numerical variables names (numeric and integer)
factor_var <- sapply(training_set, is.factor) #factor variables names


### Model building

## Full model training
# Create imputed and normalized training/test sets
#Source functions
source(paste0(dirname(dir), '/impute_mean_mode.R'))
source(paste0(dirname(dir), '/normalize_min_max.R'))
#Identify quantitative/categorical variables
categorical_var <- which(colnames(training_set) %in% c('gender','education','marital_status','smoking',
                                                        'depression_scale','mod_vig_pa','heartd_hx','hypertension','high_chol','sr_poor_health')) #indices
quantitative_var <- which(!(1:ncol(training_set) %in% c(categorical_var,23))) #indices
mode_var <- colnames(training_set)[categorical_var] #names
mean_var <- colnames(training_set)[quantitative_var] #names

#Imputation
imputed_data <- impute_mean_mode(training_set, test_set, mean_var, mode_var)
training_set_imput <- imputed_data$training_set_imput
test_set_imput <- imputed_data$test_set_imput

#Normalization
norm_var <- colnames(training_set)[numerical_var] #names of numerical variables
normalized_data <- normalize_min_max(training_set_imput, test_set_imput, norm_var)
training_set_imput_norm<- normalized_data$training_set_norm
test_set_imput_norm <- normalized_data$test_set_norm
#Normalization fail
max(training_set_imput$hdl)
max(test_set_imput$hdl)

#Training
full_model <- glm(outcome~., family='binomial', data=training_set_imput_norm)
summary(full_model)

## Feature selection
null_model <- glm(outcome ~ 1, data = training_set_imput_norm, family = "binomial") #null model (intercept only)model_range del)
model_range <- list(lower = formula(null_model), upper = formula(full_model)) #range of models to try

backward_model <- step(full_model, model_range, direction = 'backward', trace = T) #trace=T allows to trace the process
forward_model <- step(null_model, model_range, direction = 'forward', trace = F) 
stepwise_backward_model <- step(full_model, model_range, direction = 'both', trace = F) 
stepwise_forward_model <- step(null_model, model_range, direction = 'both', trace = F)
# Comparison of selected variables in each model
names(sort(backward_model$coefficients))
names(sort(forward_model$coefficients))
names(sort(stepwise_backward_model$coefficients))
names(sort(stepwise_forward_model$coefficients))

#Stability of backward selection
Nboot <- 50 #number of bootstrap iterations

features_bootstrap <- character() #initialize selected features aggregator
set.seed(123)
for(i in 1:Nboot){
  cat('Iteration number = ', i, '\n') #counter
  
  ind_0_sample <- sample(ind_0_training, length(ind_0_training), replace = T)
  ind_1_sample <- sample(ind_1_training, length(ind_1_training), replace = T)
  ind_sample <- c(ind_0_sample, ind_1_sample)
  training_set_bootstrap <- training_set[ind_sample, ]
  test_set_bootstrap <- training_set[-ind_sample, ] 
  
  #Imputation
  imputed_data <- impute_mean_mode(training_set_bootstrap, test_set_bootstrap, mean_var, mode_var)
  training_set_imput_bootstrap <- imputed_data$training_set
  test_set_imput_bootstrap <- imputed_data$test_set
  
  #Normalization
  normalized_data <- normalize_min_max(training_set_imput_bootstrap, test_set_imput_bootstrap, norm_var)
  training_set_imput_norm_bootstrap <- normalized_data$training_set_norm
  test_set_imput_norm_bootstrap <- normalized_data$test_set_norm
  
  #Define baseline models
  full_bootstrap <- glm(outcome ~ ., data = training_set_imput_norm_bootstrap, family = "binomial") #full model (all predictors)
  null_bootstrap <- glm(outcome ~ 1, data = training_set_imput_norm_bootstrap, family = "binomial") #null model (intercept only)
  range_bootstrap <- list(lower = formula(null_bootstrap), upper = formula(full_bootstrap))
  
  #Backward selection
  backward_bootstrap <- step(full_bootstrap, range_bootstrap, direction = 'backward', trace = F)
  features <- colnames(backward_bootstrap$model) #current model selected features
  features <- features[-which(features=='outcome')] #remove outcome
  
  features_bootstrap <- c(features_bootstrap, features) #concatenate in the aggregator
  #features only need to be aggregated and counted, no need to reference the model they come from
}

# Percentage of each selected feature
features_perc <- 100*table(features_bootstrap)/Nboot #features percentages
sort(features_perc, decreasing = T) #the values are coherent with percentages
# Feature selection
perc <- 80 #percentage
selected_features <- names(features_perc[features_perc>perc]) #frequency threshold
#note that feature_perc has a "names" attribute

#Train final model on the imputed normalized sets
model_formula <- as.formula(paste('outcome', paste(selected_features, collapse = '+'), sep = '~'))
#self explanatory syntax to define the model formula
final_backward_model <- glm(model_formula, training_set_imput_norm, family = 'binomial')

## LASSO
library(glmnet)
library(pROC)
K <- 5 #number of folds for cv
lambda_vec <- c(seq(0.0001,0.001,0.0001),seq(0.001,0.01,0.001),seq(0.01,0.1,0.01)) #coefficient values
auroc_lambda_average <- numeric(length(lambda_vec)) #initialize average AUROCs vector

for(i in 1:length(lambda_vec)){
  auc_cv <- numeric(length(K)) #initialize 5 fold AUCs
  
  # 5 fold CV
  ind_0_start <- ind_0_training #assign to new variable since they'll get modified
  ind_1_start <- ind_1_training #same
  
  set.seed(123)
  for (k in 1:K) {
    # Sample the indices
    ind_0_cv <- sample(ind_0_start, round(length(ind_0_training)/K), replace=F) #sampling of size #training_set/k
    ind_1_cv <- sample(ind_1_start, floor(length(ind_1_training)/K), replace=F) #use floor otherwise it exceeds the length
    ind_0_start <- setdiff(ind_0_start,ind_0_cv) #remove already extracted indices 
    ind_1_start <- setdiff(ind_1_start,ind_1_cv) #same
    
    # Create training/test sets of the current fold
    test_set_cv <- training_set[c(ind_0_cv,ind_1_cv),] 
    training_set_cv <- training_set[-c(ind_0_cv,ind_1_cv),]
    
    #Imputation
    imputed_data <- impute_mean_mode(training_set_cv, test_set_cv, mean_var, mode_var)
    training_set_imput_cv <- imputed_data$training_set
    test_set_imput_cv <- imputed_data$test_set
    
    #Normalization
    normalized_data <- normalize_min_max(training_set_imput_cv, test_set_imput_cv, norm_var)
    training_set_imput_norm_cv <- normalized_data$training_set_norm
    test_set_imput_norm_cv <- normalized_data$test_set_norm
    
    # Current fold model training
    M_training <- model.matrix(outcome ~ ., training_set_imput_norm_cv)[,-1] #matrix formula
    M_test <- model.matrix(outcome ~ ., test_set_imput_norm_cv)[,-1] #same for test
    lasso_cv <- glmnet(M_training, training_set_imput_norm_cv$outcome, family = 'binomial', 
                       alpha=1, lambda=lambda_vec[i], standardize=F)
    
    # Current fold model validation
    lasso_cv_pred <- predict(lasso_cv, newx=M_test, type='response') #predictions
    auc_cv[k] <- roc(test_set_imput_norm_cv$outcome, as.vector(lasso_cv_pred), quiet=T)$auc #AUROC
    
  }
  
  # For each lambda assign the cross validated values
  auroc_lambda_average[i] <- mean(auc_cv) #mean over 5 folds
  
}
auroc_max <- max(auroc_lambda_average)

# Optimal lambda
max_difference <- 0.005 #maximum difference between max and average of the CV
lambda_opt_ind <- which.max(lambda_vec[which(auroc_max-auroc_lambda_average<max_difference)])
lambda_opt <- lambda_vec[lambda_opt_ind] #optimal lambda with parsimony criterion
lambda_max_ind <- which.max(auroc_lambda_average)
lambda_max <- lambda_vec[lambda_max_ind] #optimal lambda as the one maximizing AUC
cat('Optimal lambda value with 5-fold Cross Validation:', lambda_opt)


# Visualize AURCOC(lambda) (zoomed in 0.0001-0.01)
plot(lambda_vec[1:19], auroc_lambda_average[1:19], type='b',
     main='Mean 5-fold CV AUROC(lambda) ', xlab='lambda', ylab='AUROC')
points(c(lambda_opt,lambda_max), c(auroc_lambda_average[lambda_opt_ind], auroc_lambda_average[lambda_max_ind]), 
       col=c('green','red'), lwd=3)
text(c(lambda_opt,lambda_max), c(auroc_lambda_average[lambda_opt_ind], auroc_lambda_average[lambda_max_ind]),
     labels=c('parsimony lambda','maximizing AUC lambda'),adj = c(0, 0))

# Train optimal model with LASSO
M_training <- model.matrix(outcome ~ ., training_set_imput_norm)[,-1] #create a matrix based formula
M_test <- model.matrix(outcome ~ ., test_set_imput_norm)[,-1] #create a matrix based formula
model_lasso_final <- glmnet(M_training, training_set_imput_norm$outcome, family = 'binomial',
                            alpha = 1, lambda = lambda_opt, standardize = F)

model_lasso_final$beta #display coefficients


### Validation
# pROC::ci() performs 2000 stratified samples bootstrap distribution estimate
par(mfrow=c(1,3))

## Full model
full_model_pred <- predict(full_model, newdata = test_set_imput_norm, type = 'response') #probabilities
roc_full_model <- roc(response = test_set_imput_norm$outcome, predictor = full_model_pred, quiet=T)
auroc_full_model <- roc_full_model$auc
ci_full_model <- ci(roc_full_model)
print(paste0('Optimal full model AUROC (95% C.I.): ', round(auroc_full_model, digits=3),
             ' (', round(ci_full_model[1], digits=3), '-', round(ci_full_model[3], digits=3), ')'))
x_full <- 1-roc_full_model$specificities #extract x axis (1-specificity)
y_full <- roc_full_model$sensitivities #extract y axis (sensitivity)
plot(x_full, y_full, main = 'Full Model', xlab = '1 - Specificity', ylab = 'Sensitivity', type = 'l') #plot
lines(c(0,1),c(0,1),col='red') #add random classifier line (diagonal)
lines(c(0,0),c(0,1), col='blue')
lines(c(1,0),c(1,1), col='blue')

# Backward model
backward_pred <- predict(backward_model, test_set_imput_norm, type = 'response')
roc_backward <- roc(test_set_imput_norm$outcome, backward_pred, quiet=T)
auroc_backward <- roc_backward$auc
ci_backward <- ci(roc_backward)
print(paste0('Optimal Backward model AUROC (95% C.I.): ', round(auroc_backward, digits=3),
             ' (', round(ci_backward[1], digits=3), '-', round(ci_backward[3], digits=3), ')'))
x_back <- 1-roc_backward$specificities #extract x axis (1-specificity)
y_back <- roc_backward$sensitivities #extract y axis (sensitivity)
plot(x_back, y_back, main = 'Backward Model', xlab = '1 - Specificity', ylab = 'Sensitivity', type = 'l') #plot
lines(c(0,1),c(0,1),col='red') #add random classifier line (diagonal)
lines(c(0,0),c(0,1), col='blue')
lines(c(1,0),c(1,1), col='blue')

# LASSO
lasso_final_pred <- predict(model_lasso_final, M_test, type = 'response')
roc_lasso <- roc(test_set_imput_norm$outcome, as.vector(lasso_final_pred), quiet=T)
auroc_lasso <- roc_lasso$auc
ci_lasso <- ci(roc_lasso)
print(paste0('Optimal LASSO model AUROC (95% C.I.): ', round(auroc_lasso, digits=3),
            ' (', round(ci_lasso[1], digits=3), '-', round(ci_lasso[3], digits=3), ')'))
x_lasso <- 1-roc_lasso$specificities #extract x axis (1-specificity)
y_lasso <- roc_lasso$sensitivities #extract y axis (sensitivity)
plot(x_lasso, y_lasso, main = 'LASSO Model', xlab = '1 - Specificity', ylab = 'Sensitivity', type = 'l') #plot
lines(c(0,1),c(0,1),col='red') #add random classifier line (diagonal)
lines(c(0,0),c(0,1), col='blue')
lines(c(1,0),c(1,1), col='blue')

#Final LASSO model predictions and accuracy
lasso_final_class <- as.factor(ifelse(lasso_final_pred>0.5,1,0))
caret::confusionMatrix(lasso_final_class, test_set_imput_norm$outcome)
acc_lasso_final <- max(table(lasso_final_class==test_set_imput_norm$outcome))/nrow(test_set_imput_norm)
acc_lasso_final <- round(acc_lasso_final, digits=4)*100

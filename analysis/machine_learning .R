library(caret)

# Read in Data 
# look at the structure of the data -- how is the dependent variable coded

# Starting with Freq 3 
str(freq.3.subset)
freq.3.ml = freq.3.subset

freq.3.ml <- within(freq.3.ml, {
  success <- factor(success)
  antigentype = NULL
  infected = NULL
  freq = NULL
  name = as.factor(name)
  netau = NULL
  simDay = NULL
})

str(freq.3.ml)
# 2. Pre-process using caret

# 2.1 check if there are missing values, you can impute using KNN
# 2.2 scale and center the numerical data 
preProcValues <- preProcess(freq.3.ml, method = c("knnImpute","center","scale"))

library(RANN)
train_processed <- predict(preProcValues, freq.3.ml)
sum(is.na(train_processed))

# now the data set I want to use is train_processed

# 2.3 Hot encode for dummy variables -- but need to change the 
#dependent variable to numerical

#train_processed$Loan_Status <- ifelse(train_processed$Loan_Status == "N", 0, 1)
train_processed$success <- ifelse(train_processed$success == "Est.", 1, 0)

#id <- train_processed$Loan_ID
#train_processed$Loan_ID <- NULL # Not sure why doing this ... will check later

# Converting every categorical variable to numerical
dmy <- dummyVars("~.", data = train_processed, fullRank = T)
train_transformed <- data.frame(predict(dmy, newdata = train_processed))

# Check the structure of the transformed file -- it should have now every categorical variable listed out
str(train_transformed)
train_processed$success=as.factor(train_processed$success)


### 3. Splitting the data set for cross-validation 
index <- createDataPartition(train_transformed$success, p = .75, list = FALSE)
trainSet <- train_transformed[index,]
testSet <- train_transformed[-index,]

trainSet$success = as.factor(trainSet$success)

### 4. Feature Selection using Caret 

# using the recursive feature to find the best subset of features
control <- rfeControl(functions = rfFuncs,
                      method = "repeatedcv",
                      repeats = 3,
                      verbose = FALSE)
outcomeName <- 'success'
#predictors <- names(trainSet)[!names(trainSet) %in% outcomeName]
predictors <- c("individual.varSigma", "meanLoad", "dominant.freq",  "diversity", "tmrca", "serialInterval",
                "antigenicDiversity", "totalS", "meanR", "varR", "varBeta", "varSigma", "entropy", "ratio.mutation", "ratio.meanR",
                "ratio.meanBeta", "ratio.varBeta", "ratio.varSigma")

success_Pred_Profile <- rfe(trainSet[,predictors], trainSet[,outcomeName],
                         rfeControl = control)
success_Pred_Profile

# taking only the top 5 
#predictors <- c("Credit_History", "LoanAmount", "Loan_Amount_Term", "ApplicantIncome", "CoapplicantIncome")


# 5. Training models using Caret

fitControl <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats =5)


model_gbm<-train(trainSet[,predictors],trainSet[,outcomeName],method='gbm',
                 trControl = fitControl)
#model_rf<-train(trainSet[,predictors],trainSet[,outcomeName],method='rf')
model_nnet<-train(trainSet[,predictors],trainSet[,outcomeName],method='nnet', trControl = fitControl)
model_glm<-train(trainSet[,predictors],trainSet[,outcomeName],method='glm', trControl = fitControl)

#6. Parameter tuning 
# 5-fold cross-validation repreated 5 times


#7 Variable importance estimation using caret
varImp(object = model_glm)
varImp(object = model_gbm)
varImp(object = model_nnet)

plot(varImp(object=model_gbm),main="GBM - Variable Importance")
# Can use this technique to check the importance of different models


#8. Predictions using Caret

predictions <- predict.train(object = model_glm, testSet[,predictors], type = "raw")
table(predictions)
confusionMatrix(predictions,testSet[,outcomeName])

#9. Redoing with a smaller subset of samples 

predictors = c("ratio.meanR", "varR", "entropy", "ratio.mutation", "antigenicDiversity")
model_gbm<-train(trainSet[,predictors],trainSet[,outcomeName],method='gbm',trControl = fitControl)
model_rf<-train(trainSet[,predictors],trainSet[,outcomeName],method='rf',  trControl = fitControl)
model_nnet<-train(trainSet[,predictors],trainSet[,outcomeName],method='nnet', trControl = fitControl)
model_glm<-train(trainSet[,predictors],trainSet[,outcomeName],method='glm',  trControl = fitControl)

predictions <- predict.train(object = model_gbm, testSet[,predictors], type = "raw")
confusionMatrix(predictions,testSet[,outcomeName])



##### Feature Selection Tutorial
library('Metrics')
library('randomForest')
library('ggthemes')
library(dplyr)

set.seed(101)
data <- read.csv("~/Downloads/stock_data.csv", stringsAsFactors = T)
data$Y = as.factor(data$Y) # specifcying outcome variable as factor
data$Time <- NULL

# Dividing the dataset into train and test

train <- data[1:2000,]
test <- data[2001:3000,]

# applying Random Forest

model_rf <-randomForest(Y~., data = train)
preds <- predict(model_rf, test[,-101])
table(preds)

# checking accuracy
auc(preds, test$Y)

variable.import = data.frame(importance(model_rf))

# Pick a top choice

preds <- predict(model_rf, test[,-101])
table(preds)
auc(preds, test$Y)

plot(varImp(object=model_rf),main="RF - Variable Importance")

varImp(object = model_rf)

########################### Trying with a Decision Tree 
library(caret)
library(rpart.plot)

dim(testSet); dim(trainSet)

model_tree <-train(trainSet[,predictors],trainSet[,outcomeName],method='rpart',
                  trControl = fitControl, parms = list(split = "information"),
                  tuneLength = 10)
model_tree

prp(model_tree$finalModel, box.palette = "Reds", tweak = 1.2)

predictions <- predict.train(object = model_tree, testSet[,predictors], type = "raw")
confusionMatrix(predictions, testSet$success, positive = "1")


# collect resamples
results <- resamples(list(GBM = model_gbm, GLM=model_glm, 
                          NNET=model_nnet, RF=model_rf,
                          DT=model_tree))
# summarize the distributions
summary(results)
# boxplots of results
bwplot(results)
# dot plots of results
dotplot(results)

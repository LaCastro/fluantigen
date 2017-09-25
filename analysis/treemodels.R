excluded = c("postAntigen", "day", "simDay", "date", "totalN", "totalR",
             "totalCases", "dominant.type", "day.1")
predictorNames = colnames(freq.two)


frequencysub = freq.five[, -which(variableNames%in%excluded)]
frequencysub$success2 = NA
frequencysub$success2[which(frequencysub$success == "yes")] = 0
frequencysub$success2[which(frequencysub$success == "no")] = 1

str(frequencysub)
# Creating Dummy variables 
frequencysub$.id = as.factor(frequencysub$.id)
frequencysubDummy <- dummyVars("~.", data = frequencysub, fullRank = F)
frequencysub <- as.data.frame(predict(frequencysubDummy, frequencysub))

# if the proportion was smaller than 15% would be considered a rare event and would be more challening to model
# that's why check it at the 5/10%
prop.table(table(frequencysub$success2))

outcomeName <- 'success'
predictorsNames <- names(frequencysub)[names(frequencysub)!=outcomeName]
predictorsNames=predictorsNames[-21]

# For GBM - classification 
frequencysub$success2 <- ifelse(frequencysub$success==2,'yes','nope')
frequencysub$success <- as.factor(frequencysub$success)
outcomeName <- 'success'

set.seed(1234)
splitIndex <- createDataPartition(frequencysub[,outcomeName], p = .75, list = FALSE, times = 1)
trainDF <- frequencysub[ splitIndex,]
testDF  <- frequencysub[-splitIndex,]

# Control the resampling of your data - splitting the training set internally to figure out the best settings for the model
objControl <- trainControl(method='cv', number=3, returnResamp='none', summaryFunction = twoClassSummary, classProbs = TRUE)

objModel <- train(trainDF[,predictorsNames], trainDF[,outcomeName], 
                  method='gbm', 
                  trControl=objControl,  
                  metric = "ROC",
                  preProc = c("center", "scale"))

summary(objModel)

# Class Predictions
predictions <- predict(object = objModel, testDF[,predictorsNames], type = 'raw')
print(postResample(pred=predictions, obs=as.factor(testDF[,outcomeName])))


# Probabilites 
predictions <- predict(object=objModel, testDF[,predictorsNames], type='prob')
auc <- roc(ifelse(testDF[,outcomeName]=="yes",1,0), predictions[[2]])
print(auc$auc)


#### Regression
outcomeName <- 'success'
set.seed(1234)
splitIndex <- createDataPartition(frequencysub[,outcomeName], p = .8, list = FALSE, times = 1)
trainDF <- frequencysubsub[ splitIndex,]
testDF  <- frequencysubsub[-splitIndex,]

objControl <- trainControl(method='cv', number=3, returnResamp='none')
objModel <- train(trainDF[,predictorsNames], trainDF[,outcomeName], method='glmnet',  metric = "RMSE", trControl=objControl)

predictions <- predict(object=objModel, testDF[,predictorsNames])
auc <- roc(testDF[,outcomeName], predictions)
print(auc$auc)

plot(varImp(objModel, scale = F))



################ Logistic Regression
library(pscl)
library(ROCR)
library(caret)

str(freq.five)
freq.five$.id = as.factor(freq.five$.id)
freq.five$totalI=as.numeric(freq.five$totalI)

str(freq.two)

variableNames = colnames(freq.five)
excluded = c("postAntigen", "day", "simDay", "date", "totalN", "totalR", "totalI", 
             "totalCases", "dominant.type", "day.1", "max.I", "min.I")

frequencysub = freq.five[, -which(variableNames%in%excluded)]
prop.table(table(frequencysub$success))
outcomeName = "success"

set.seed(505)
splitIndex <- createDataPartition(frequencysub[,outcomeName], p = .80, list = FALSE, times = 1)
train <- frequencysub[ splitIndex,]
test  <- frequencysub[-splitIndex,]

### Model Fitting 
model <- glm(success~., family = binomial(link = 'logit'), data = train)
summary(model)
anova(model, test="Chisq")

drop1(model, test="Chisq")
backwards = step(model, direction="backward")

add1(backwards, ~.^2,test="Chisq")
interaction = step(backwards, ~.^2, direction = "forward")

search$anova

# Model Fit 
pR2(search)

test$success = as.factor(as.numeric(test$success))

fittedResults <- predict(interaction, newdata = test[,-1], type = 'response')
fittedResults <- ifelse(fittedResults > 0.5,1,0)

test$success = ifelse(test$success == "yes",1,0)

misClasificError <- mean(fittedResults != test[,1])
print(paste('Accuracy', 1-misClasificError))

p <- predict(interaction, newdata = subset(test[,-1]), type = "response")
pr <- prediction(p, test$success)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)

auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc

step(model, ~.^2, direction = "both")
anova(backwards, test = "Chisq")

install.packages('lme4')

glmer()

library(caret)
library(randomForest)

table(dataPurged$success, dataPurged$quarter)

dataPurged$transition = NA
dataPurged[dataPurged$success==1, "transition"] = "yes"
dataPurged[dataPurged$success==0, "transition"] = "no"
dataPurged$transition=as.factor(dataPurged$transition)

featurePlot(x=dataPurged[,3:10],
            y=dataPurged$transition,
            plot="box",
            scales = list(x=list(relation="free"),
                          y=list(relation="free")))


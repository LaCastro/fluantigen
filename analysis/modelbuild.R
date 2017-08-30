library(caret)
library(pROC)


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




###########
hdp <- read.csv("https://stats.idre.ucla.edu/stat/data/hdp.csv")

hdp <- within(hdp, {
  Married <- factor(Married, levels = 0:1, labels = c("no", "yes"))
  DID <- factor(DID)
  HID <- factor(HID)
})

m  <- glmer(remission~IL6+CRP+CancerStage+LengthofStay+Experience+(1|DID), data = hdp, family = binomial,
            control = glmerControl(optimizer = "bobyqa"), nAGQ=10)

summary(m)
print(m, corr=FALSE)
library(car)
Anova(m)
head(hdp)

#########################################################################################
predictorNames = colnames(freq.five)

#excluded = c("postAntigen", "day", "simDay", "date", "totalN", "totalR",
#             "totalCases", "dominant.type", "max.I", "min.I")

included = c("diversity", "meanBeta", "covBetaSigma", "dominant.freq", "normalize.I", "success", ".id")

data = freq.five[, which(predictorNames%in%included)]

data <- within(data, {
  success <- factor(success, levels=c("yes", "no"), labels = c(1,0))
  .id <- factor(.id)
})

set.seed(1234)
outcomeName = "success"
splitIndex <- createDataPartition(data[,outcomeName], p = .75, list = FALSE, times = 1)
trainDF <- data[ splitIndex,]
testDF  <- data[-splitIndex,]


pvars <- included[1:5]

dataS <- trainDF
dataS[pvars] <- lapply(dataS[pvars],scale)

testDF[pvars] <- lapply(testDF[pvars],scale)


trial.id <- glmer(success~1 +(1|.id), data = dataS, family = binomial,
                  control = glmerControl(optimizer="bobyqa"),nAGQ=10)


trial.meanBeta <- glmer(success~meanBeta +(1|.id), data = dataS, family = binomial,
                        control = glmerControl(optimizer="bobyqa"),nAGQ=10)

anova(trial.id, trial.meanBeta)

trial.meanBetafreq <- glmer(success~meanBeta+dominant.freq + (1|.id), data = dataS, 
                            family = binomial, control = glmerControl(optimizer="bobyqa"), nAGQ=10)

trial.meanBetaNormalizeI <- glmer(success~meanBeta+normalize.I + (1|.id), data = dataS, 
                                  family = binomial, control = glmerControl(optimizer="bobyqa"), nAGQ=10)

trial.full <- glmer(success~meanBeta+normalize.I+dominant.freq + (1 |.id), data = dataS,
                    family = binomial, control = glmerControl(optimizer="bobyqa"), nAG=10)
summary(trial.full)


exp(predict(trial.full, newdata = testDF))

trainDF %>%
  gather(key = variable, value = value, -success,-.id) %>%
  ggplot(aes(value))+geom_density(aes(color = success)) +
  facet_wrap(~variable,scales = "free")


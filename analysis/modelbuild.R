library(caret)
library(pROC)
library(lme4)


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

excluded = c("postAntigen", "day.1", "simDay", "date", "totalN", "totalR",
             "totalCases", "dominant.type", "totalI", "max.I", "min.I")

#included = c("diversity", "meanBeta", "covBetaSigma", "dominant.freq", "normalize.I", "success", ".id")

data = freq.five[, -which(predictorNames%in%excluded)]

str(data)
data <- within(data, {
  success <- factor(success, levels=c("yes", "no"), labels = c(1,0))
  .id <- factor(.id)
  day <- as.numeric(as.character(day))
})

set.seed(1234)

##### 
outcomeName = "success"

# Will do the splitting when I have more trials 
splitIndex <- createDataPartition(data[,outcomeName], p = .75, list = FALSE, times = 1)
dataScaled <- data
dataScaled[,-c(1, 21)] <- lapply(dataScaled[,-c(1,21)],scale)

trainDF <- dataScaled[ splitIndex,]
testDF  <- dataScaled[-splitIndex,]



### Testing for multicollinearity 
trial.all <- glmer(success~.-.id + (1|.id), data = trainDF, family = binomial,
                        control = glmerControl(optimizer="bobyqa"),nAGQ=10)
vif.mer(trial.all)

####### Based on variables that are not colinear
modelpurged <- glmer(success ~ diversity + dominant.freq+antigenicTypes+serialInterval+
                       antigenicDiversity+netau+normalize.I+meanLoad+tmrca + (1|.id),
                     data = dataScaled, family = binomial,
                     control = glmerControl(optimizer="bobyqa"),nAGQ=10)


## going to exclude the data of the other variables 
vif.included = c("success", ".id", "diversity", "dominant.freq", "antigenicTypes", "serialInterval", "antigenicDiversity", "netau", "normalize.I", "meanLoad", "tmrca")
trainDF = dataScaled[splitIndex, which(colnames(dataScaled) %in% vif.included)]

### Null Model
modelnull = glmer(success ~ 1 + (1|.id),
                  data = trainDF, family = binomial,
                  control = glmerControl(optimizer="bobyqa"),nAGQ=10)
summary(modelnull)



# Building Models with one term 
model1 <- glmer(success ~ meanLoad + (1|.id),
                data = trainDF, family = binomial,
                control = glmerControl(optimizer="bobyqa"),nAGQ=10)

summary(model1)
anova(modelnull,model1)

model2 <- glmer(success ~ normalize.I + (1|.id),
                     data = trainDF, family = binomial,
                     control = glmerControl(optimizer="bobyqa"),nAGQ=10)
summary(model2)
anova(model1, model2)

# Building models with 2 terms -- null model becomes the one with meanLoad
model1 <- glmer(success ~ meanLoad + dominant.freq + (1|.id),
                data = trainDF, family = binomial,
                control = glmerControl(optimizer="bobyqa"),nAGQ=10)
anova(modelnull, model1)
summary(model1)

model2 <- glmer(success ~ meanLoad + normalize.I + (1|.id),
                data = trainDF, family = binomial,
                control = glmerControl(optimizer="bobyqa"),nAGQ=10)
summary(model2)
anova(model1, model2)

############# Building model 3 terms -- null model becomes one with meanLoad and Dolminant Frequency
model1 <- glmer(success ~ meanLoad + antigenicTypes + dominant.freq + (1|.id),
                data = trainDF, family = binomial,
                control = glmerControl(optimizer="bobyqa"),nAGQ=10)
summary(model1)
anova(modelnull, model1)

model2 <- glmer(success ~ meanLoad + antigenicTypes + normalize.I + (1|.id),
                data = trainDF, family = binomial,
                control = glmerControl(optimizer="bobyqa"),nAGQ=10)
summary(model2)

################ Building model 4 terms -- null model becomes meanLoad/dominantfreq/normalize.I + (1|.id)
model1 <- glmer(success ~ meanLoad + dominant.freq + normalize.I + antigenicTypes + (1|.id),
                data = trainDF, family = binomial,
                control = glmerControl(optimizer="bobyqa"),nAGQ=10)
summary(model1)
anova(modelnull, model1)


############### looking at interaction terms
model1 <- glmer(success ~ (meanLoad + dominant.freq + antigenicTypes)^2 + (1|.id),
                data = trainDF, family = binomial,
                control = glmerControl(optimizer="bobyqa"),nAGQ=10)
summary(model1)
anova(modelnull, model1)


               
############### testing prediction
library(arm)
library(pROC)
model.predict = predict(modelnull, newdata = testDF, type = "response")
glmer.ROC <- roc(predictor=model.predict, response=testDF$success)
glmer.ROC$auc

plot(glmer.ROC)


plot(dataScaled$meanLoad, dataScaled$success, ylab = "Probability of Success")
binnedplot(fitted(modelnull),resid(modelnull))


############### Without Tropics effects 
modelbase <- glm(success~., family = binomial(link = 'logit'), data = dataVif)
summary(modelbase)
anova(modelbase, test="Chisq")

lme.bestfit = step(modelbase,  ~., test="Chisq",direction = "backward", data = dataVif)

fittedResults = predict(lme.bestfit, newdata = dataVif, type = "response")
lme.roc <- roc(predictor=fittedResults, response=dataVif$success, )
lme.roc


#############
library(sjPlot)
library(sjmisc)

data
set_theme(theme = "forest", 
          geom.label.size = 3, 
          axis.textsize = .9, 
          axis.title.size = .9)

# create binary response
efc$hi_qol <- dicho(efc$quol_5)
# prepare group variable
efc$grp = as.factor(efc$e15relat)
levels(x = efc$grp) <- get_labels(efc$e15relat)
# data frame for fitted model
mydf <- data.frame(hi_qol = efc$hi_qol,
                   sex = to_factor(efc$c161sex),
                   c12hour = efc$c12hour,
                   neg_c_7 = efc$neg_c_7,
                   grp = efc$grp)
# fit glmer
fit2 <- glmer(hi_qol ~ sex + c12hour + neg_c_7 + (1 | grp),
              data = mydf, family = binomial("logit"))

sjp.glmer(modelnull, type = "fe")
plot(modelnull)
library(ggplot2)

ggplot(data.frame(eta=predict(modelnull,type="link"),pearson=residuals(modelnull,type="pearson")),
       aes(x=eta,y=pearson)) +
  geom_point() +
  theme_bw()

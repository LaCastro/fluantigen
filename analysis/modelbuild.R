library(plyr)
library(dplyr)
library(tidyverse)
library(lme4)
library(pROC)
library(caret)
library(cvTools)

#########################################################################################
predictorNames = colnames(antigen.data)

exclude.emerge = c("oriAntigen", "postAntigen", "cases", "cumulativeTypes", 
                   "final.max", "dominant.type", "life.length", "max.I", "min.I",
                   "days.above", "infected", "N", "R")

data = antigen.data[, -which(predictorNames%in%exclude.emerge)]


predictorNames = colnames(freq.05)
excluded = c("postAntigen", "oriAntigen",  "day.1", "date", "totalN", "totalR",
             "totalCases", "dominant.type", "totalI", "max.I", "min.I")
data = freq.05[, -which(predictorNames%in%excluded)]

data <- within(data, {
  success <- factor(success, levels=c("yes", "no"), labels = c(1,0))
  .id <- factor(.id)
  day <- as.numeric(as.character(day))
  #mutLoad <- as.numeric(as.character(mutLoad))
  #distance <- as.numeric(as.character(distance))
  simDay <- as.numeric(as.character(simDay))
})

#data %>%
#  filter(day > (3*365) & day < 3650) -> dataBurn

data %>%
  filter(day > 365) -> dataBurn

number.years.rep = 17

##### When the simDay isn't in terms of the burn in day, don't need the -1 and -3 
dataBurn %>%
  mutate(time.of.year = (simDay)-floor((x = day/365))) %>%
  mutate(quarter = ifelse(time.of.year < .25, 1,
                          ifelse(time.of.year < .5, 2,
                                 ifelse(time.of.year < .75, 3,4)))) -> dataScaled

dataScaled$quarter = as.factor(dataScaled$quarter)

# Take our any data that's before 3 years
table(dataScaled$success)/nrow(dataScaled)

dataScaled%>%
  group_by(.id)%>%
  summarize(mutations = n())%>%
  summarize(mean.mutations = mean(mutations)/number.years.rep)


factor.variables = which(sapply(dataScaled,is.factor)==TRUE)
dataScaled$netau[!is.finite(dataScaled$netau)] <- NA

dataScaled[,-factor.variables] <- lapply(dataScaled[,-factor.variables],scale)


trial.all <- glmer(success ~.-.id -quarter+ (1|.id), data = dataScaled, family = binomial,
                   control = glmerControl(optimizer="bobyqa"),nAGQ=0)

vif.mer(trial.all)

vif.included = c("success",
                 "quarter",
                 "distance",
                 "mutLoad",
                 "diversity", 
                 "tmrca",
                 "netau",
                 "serialInterval",
                 "antigenicDiversity",
                 "dominant.freq",
                 "meanLoad",
                 "antigenicTypes",
                 "normalize.I",
                 "time.of.year",  # 
                 "varR",
                 "meanR",
                 "covBetaSigma")

dataPurged = dataScaled[,which(colnames(dataScaled) %in% vif.included)]
dataPurged = cbind(dataScaled[,".id"], dataPurged)
colnames(dataPurged)[1] = "trial"

trial.purged <- lme4::glmer(success ~.-trial -quarter+ (1 | trial), data = dataPurged, family = binomial,
                   control = glmerControl(optimizer="bobyqa"),nAGQ=0)

vif.mer(trial.purged)

################### Doing K-cross validation 
outcomeName = "success"

set.seed(62417)
folds.length = 5;
folds = cvFolds(n = length(unique(dataPurged$trial)), K = folds.length, type = "random")

glmerperf=rep(0, folds.length); glmperf=glmerperf;
aic = rep(0,folds.length)
coef = rep(0,folds.length)
ref.sd = rep(0,folds.length)

dataPurged %>%
  gather(key = variable, value = value, -trial, -success) -> data.l

data.l$value = as.numeric(as.character(data.l$value))
str(data.l)

single.model.results = ddply(data.l, .var = c("variable"), .fun = function(x) {
  # K-fold validation for each model based on a single term 
  for(n in 1:folds.length) { 
    test.ids = folds$subsets[folds$which==n]
    test.trials = unique(data.l$trial)[test.ids]
    
    testdata = x[which(x$trial %in% test.trials), ]
    traindata = x[-which(x$trial %in% test.trials),]
    
    GLMER <- lme4::glmer(success ~ value  + (1 | trial), data = traindata, family="binomial", 
                         control = glmerControl(optimizer="bobyqa"),nAGQ=0)
   
    aic[n] <- extractAIC(GLMER)[2]
    coef[n] <- summary(GLMER)$coefficients[2,4]
    ref.sd[n] <- sqrt(VarCorr(GLMER)$trial[1])
    glmer.probs <- predict(GLMER, newdata=testdata, type="response", allow.new.levels=TRUE)
    glmer.ROC <- roc(predictor=glmer.probs, response=testdata$success)
    glmerperf[n] <- glmer.ROC$auc
  }
  single.performance=mean(glmerperf)
  single.aic=mean(aic)
#  single.coef=mean(coef)
  single.ref.sd=mean(ref.sd)
  print(x$variable[1])
  return(cbind(single.performance, single.aic, single.ref.sd))
})

single.model.results


################################## Running cross validation on single model
 
set.seed(62417)
folds.length = 5;
folds = cvFolds(n = length(unique(dataPurged$trial)), K = folds.length, type = "random")

glmerperf=rep(0, folds.length); glmperf=glmerperf;
aic = rep(0,folds.length)



for(n in 1:folds.length) { 
  test.ids = folds$subsets[folds$which==n]
  test.trials = unique(dataPurged$trial)[test.ids]
  
  testdata = dataPurged[which(dataPurged$trial %in% test.trials), ]
  traindata = dataPurged[-which(dataPurged$trial %in% test.trials),]
  
  
#  model.emerge <- glm(success~distance+mutLoad + normalize.I + meanLoad+ meanLoad*mutLoad + mutLoad*normalize.I,
#                           family = "binomial", data = traindata)
  
# model.two <- glm(success~normalize.I +covBetaSigma+serialInterval, 
#                   family = "binomial", data = traindata)
  
#  model.five <- glm(success~netau+dominant.freq+normalize.I,
#                   family = "binomial", data = traindata)

  model.ten <- glm(success~normalize.I+dominant.freq+varR+normalize.I*varR,
                   family = "binomial",data=traindata)  
  #model.glm <- lme4::glmer(success~distance+mutLoad+meanLoad+normalize.I + mutLoad*normalize.I + mutLoad*meanLoad + (1|trial),
  #                  family = "binomial", data = traindata, control = glmerControl(optimizer="bobyqa"))
  
  aic[n] <- extractAIC(model.ten)[2]
  #glmer.probs <- predict(model.glm, newdata=testdata, type="response", allow.new.levels=TRUE)
  glmer.probs <- predict(model.ten, newdata=testdata, type="response")
  
  glmer.ROC <- roc(predictor=glmer.probs, response=testdata$success)
  glmerperf[n] <- glmer.ROC$auc
}



########################### second.model.results 
dataPurged %>%
  gather(key = variable, value = value, -trial, -success, -netau) -> data.l
data.l$value=as.numeric(data.l$value)


double.model.results = ddply(data.l, .var = c("variable"), .fun = function(x) {
  # K-fold validation for each model based on a single term 
  for(n in 1:folds.length) { 
    test.ids = folds$subsets[folds$which==n]
    test.trials = unique(data.l$trial)[test.ids]
    
    testdata = x[which(x$trial %in% test.trials), ]
    traindata = x[-which(x$trial %in% test.trials),]
    
    GLMER <- lme4::glmer(success ~ netau + value  + (1 | trial), data = traindata, family="binomial", 
                         control = glmerControl(optimizer="bobyqa"),nAGQ=10)
    
    #GLMER <- glm(success ~ normalize.I + value, data = traindata, family="binomial")
    aic[n] <- extractAIC(GLMER)[2]
    coef[n] <- summary(GLMER)$coefficients[3,4]
    ref.sd[n] <- sqrt(VarCorr(GLMER)$trial[1])
    glmer.probs <- predict(GLMER, newdata=testdata, type="response", allow.new.levels=TRUE)
    glmer.ROC <- roc(predictor=glmer.probs, response=testdata$success)
    glmerperf[n] <- glmer.ROC$auc
  }
  double.performance=mean(glmerperf)
  double.aic=mean(aic)
  double.coef=mean(coef)
  double.ref.sd=mean(ref.sd)
  print(x$variable[1])
  return(cbind(double.performance, double.aic, double.coef, double.ref.sd))
})

double.model.results


############################### third term results
dataPurged %>%
  gather(key = variable, value = value, -trial, -success, -dominant.freq, -netau) -> data.l
data.l$value=as.numeric(data.l$value)

third.model.results = ddply(data.l, .var = c("variable"), .fun = function(x) {
  # K-fold validation for each model based on a single term 
  for(n in 1:folds.length) { 
    test.ids = folds$subsets[folds$which==n]
    test.trials = unique(data.l$trial)[test.ids]
    
    testdata = x[which(x$trial %in% test.trials), ]
    traindata = x[-which(x$trial %in% test.trials),]
    
    GLMER <- lme4::glmer(success ~ dominant.freq + netau+ value  + (1 | trial), data = traindata, family="binomial", 
                         control = glmerControl(optimizer="bobyqa"),nAGQ=10)
    #GLMER <- glm(success ~ normalize.I + covBetaSigma + value, data = traindata, family="binomial")
    
    aic[n] <- extractAIC(GLMER)[2]
    coef[n] <- summary(GLMER)$coefficients[4,4]
    ref.sd[n] <- sqrt(VarCorr(GLMER)$trial[1])
    glmer.probs <- predict(GLMER, newdata=testdata, type="response", allow.new.levels=TRUE)
    glmer.ROC <- roc(predictor=glmer.probs, response=testdata$success)
    glmerperf[n] <- glmer.ROC$auc
  }
  performance=mean(glmerperf)
  aic=mean(aic)
  coeff.sig = mean(coef)
  rd.eff=mean(ref.sd)
  print(x$variable[1])
  return(cbind(performance, aic, coeff.sig, rd.eff))
})


third.model.results

############################### fourth term results
dataPurged %>%
  gather(key = variable, value = value, -trial, -success, 
         -dominant.freq, -netau, -varR, -normalize.I,-tmrca) -> data.l
data.l$value=as.numeric(data.l$value)
# Remember to change the coefficien

six.model.results = ddply(data.l, .var = c("variable"), .fun = function(x) {
  # K-fold validation for each model based on a single term 
  for(n in 1:folds.length) { 
    test.ids = folds$subsets[folds$which==n]
    test.trials = unique(data.l$trial)[test.ids]
    
    testdata = x[which(x$trial %in% test.trials), ]
    traindata = x[-which(x$trial %in% test.trials),]
    
     GLMER <- lme4::glmer(success ~   dominant.freq + netau + varR + normalize.I+value +tmrca+(1 | trial), data = traindata, family="binomial", 
                           control = glmerControl(optimizer="bobyqa"),nAGQ=10)
    #GLMER <- glm(success ~ diversity + antigenicTypes + covBetaSigma + value, data = traindata, family="binomial")
    
    aic[n] <- extractAIC(GLMER)[2]
    coef[n] <- summary(GLMER)$coefficients[7,4]
    ref.sd[n] <- sqrt(VarCorr(GLMER)$trial[1])
    glmer.probs <- predict(GLMER, newdata=testdata, type="response", allow.new.levels=TRUE)
    glmer.ROC <- roc(predictor=glmer.probs, response=testdata$success)
    glmerperf[n] <- glmer.ROC$auc
  }
  performance=mean(glmerperf)
  aic=mean(aic)
  coeff.sig = mean(coef)
  rd.eff=mean(ref.sd)
  print(x$variable[1])
  return(cbind(performance, aic, coeff.sig, rd.eff))
  return(cbind(performance, aic))
})

six.model.results

#################### anova tests
dataPurged %>%
  select(trial,success, netau, dominant.freq, varR, normalize.I, tmrca, meanR) -> data.l

for(n in 1:folds.length) { 
  test.ids = folds$subsets[folds$which==n]
  test.trials = unique(data.l$trial)[test.ids]
  
  #  testdata = x[which(x$trial %in% test.trials), ]
  #  traindata = x[-which(x$trial %in% test.trials),]
  
  testdata = data.l[which(data.l$trial %in% test.trials), ]
  traindata = data.l[-which(data.l$trial %in% test.trials),]
  
  model.null <- lme4::glmer(success ~ netau+dominant.freq + varR+normalize.I+  (1 | trial),
                                                   data = traindata, family="binomial", 
                                                   control = glmerControl(optimizer="bobyqa"),nAGQ=10)
  
  model.test <- lme4::glmer(success ~ netau+dominant.freq + varR+normalize.I + tmrca  + (1 | trial),
                            data = traindata, family="binomial", 
                            control = glmerControl(optimizer="bobyqa"),nAGQ=10)
  
  print(anova(model.null, model.test, test = "Chisq"))
}




######################### Interaction Results 
dataPurged %>%
  #gather(key = variable, value = value, -trial, -success, -normalize.I, -covBetaSigma, antigenicTypes)  -> data.l
  select(trial,success, netau, dominant.freq, varR, normalize.I, tmrca) -> data.l

data.l$value=as.numeric(data.l$value)

coeff=matrix(nrow = folds.length, ncol = 3)
coeff=rep(NA,folds.length)

for(n in 1:folds.length) { 
    test.ids = folds$subsets[folds$which==n]
    test.trials = unique(data.l$trial)[test.ids]
    
    testdata = data.l[which(data.l$trial %in% test.trials), ]
    traindata = data.l[-which(data.l$trial %in% test.trials),]

    model.null <- lme4::glmer(success ~ netau+dominant.freq + varR+normalize.I + tmrca  + 
                                normalize.I*varR + (1 | trial),
                          data = traindata, family="binomial", 
                          control = glmerControl(optimizer="bobyqa"),nAGQ=10)
  
    print(summary(model.null))
    model.int <- lme4::glmer(success ~ netau+dominant.freq + varR+normalize.I + tmrca + 
                               normalize.I*varR + varR*tmrca + (1 | trial),
                              data = traindata, family="binomial", 
                              control = glmerControl(optimizer="bobyqa"),nAGQ=10)
    
  # coeff[n,] <- summary(model.9)[[10]][9:11,4]
    
    
    coef[n] <- summary(model.int)$coefficients[6,4]
    ref.sd[n] <- sqrt(VarCorr(model.int)$trial[1])
    aic[n] <- extractAIC(model.int)[2]
    glmer.probs <- predict(model.int, newdata=testdata, type="response", allow.new.levels=TRUE)
    glmer.ROC <- roc(predictor=glmer.probs, response=testdata$success)
    glmerperf[n] <- glmer.ROC$auc

    #print(anova(model.null, model.int, test="Chisq"))
  }

mean(glmerperf)
mean(aic)
mean(coef)
mean(ref.sd)
#colMeans(coeff)
               

plot(glmer.ROC)
############### testing prediction
library(arm)
library(pROC)
model.predict = predict(modelnull, newdata = testDF, type = "response")
glmer.ROC <- roc(predictor=model.predict, response=testDF$success)
glmer.ROC$auc

plot(glmer.ROC)
cbind(glmer.ROC$sensitivities, glmer.ROC$specificities, glmer.ROC$thresholds)





plot(dataScaled$meanLoad, dataScaled$success, ylab = "Probability of Success")
binnedplot(fitted(modelnull),resid(modelnull))



###### Visualization  
GLMER <- lme4::glmer(success ~ antigenicTypes + diversity  + (1 | trial), data = dataPurged, family="binomial", 
                     control = glmerControl(optimizer="bobyqa"),nAGQ=10)


library(sjPlot)
library(sjmisc)
library(sjlabelled)

sjp.glmer(model.glm)
lme4::fixef(GLMER)

summary(GLMER)

sjp.glmer(GLMER, 
          sort.est = "sort.all",
          type = "pred",
          vars = "dominant.freq") 
y.offset = .4)
          #vars = "antigenicTypes",
          #type = "fe.slope")

sjp.glmer(GLMER, type = "fe.slope")
sjp.int(GLMER, type = "cond")




  
  ##################### Random Forests Model ##########################
library(rpart)
library(e1071)
library(rpart.plot)
library(caret)

fitControl <- trainControl(method = "cv", number = 5)

cartGrid <- expand.grid(.cp=(1:50)*0.01)

#decision tree
tree_model <- train(success ~., data = traindata, method = "rpart", trControl = fitControl, tuneGrid = cartGrid)
print(tree_model)


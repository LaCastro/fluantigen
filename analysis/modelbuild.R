library(plyr)
library(dplyr)
library(tidyverse)
library(lme4)
library(pROC)
library(caret)
library(cvTools)
library(cowplot)

#########################################################################################
predictorNames = colnames(freq.first.check)
exclude.emerge = c("freq", "day", "cases", "dominant.type", "cumulativeTypes", "antigentype","N",
                   "totalS", "totalI", "I", "totalCases", "simDay")
                   
data = freq.first.check[, -which(predictorNames%in%exclude.emerge)]

data$success = as.factor(data$success)
data$name = as.factor(data$name)

data <- within(data, {
  success <- factor(success, levels=c("Est.", "Transient"), labels = c(1,0))
  name <- as.factor(name)
  #day <- as.numeric(as.character(day))
#  selected.antigen <- as.factor(selected.antigen)
#  mutLoad <- as.numeric(as.character(mutLoad))
#  distance <- as.numeric(as.character(distance))
#  simDay <- as.numeric(as.character(simDay))
})

# Take our any data that's before 3 years
table(data$success)/nrow(data)
dataScaled = data

dataScaled%>%
  group_by(.id)%>%
  summarize(mutations = n())%>%
  summarize(mean.mutations = mean(mutations)/number.years.rep)


factor.variables = which(sapply(dataScaled,is.factor)==TRUE)
dataScaled$netau[!is.finite(dataScaled$netau)] <- NA
dataScaled[,-factor.variables] <- lapply(dataScaled[,-factor.variables],scale)


trial.all <- lme4::glmer(success ~.-name + (1|name), data = dataScaled, family = binomial,
                   control = glmerControl(optimizer="bobyqa"),nAGQ=0)
summary(trial.all)
vif.results = data.frame(vif.mer(trial.all))

vif.excluded = c("success",
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

vif.excluded = c("first_totalS_diff", "first_totalS_ratio", "first_totalS_sd", 
                 "second_totalS_diff", "second_totalS_ratio", "second_totalS_sd",
                 "first_meanSigma_diff", "first_meanSigma_ratio", "first_meanSigma_sd",
                 "second_meanSigma_diff", "second_meanSigma_ratio", "second_meanSigma_sd",
                 "first_varSigma_diff", "first_varSigma_ratio", "first_varSigma_sd",
                 "second_varSigma_diff", "second_varSigma_ratio", "second_varSigma_sd",
                 "first_meanR_diff", "first_meanR_ratio", "first_meanR_sd",
                 "second_meanR_diff", "second_meanR_ratio", "second_meanR_sd")

#dataPurged = dataScaled[,which(colnames(dataScaled) != "meanSigma")]
dataPurged = data.wide[,-which(colnames(data.wide) %in% vif.excluded)]
colnames(dataPurged)[1] = "trial"

trial.purged <- lme4::glmer(success ~.-trial-selected.antigen  + (1 | trial), data = dataPurged, family = binomial,
                   control = glmerControl(optimizer="bobyqa"),nAGQ=0)
vif.results = data.frame(vif.mer(trial.purged))

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
  ungroup() -> dataPurged
dataPurged%>%
  gather(key = variable, value = value, 4:ncol(dataPurged)) -> data.scaled.l


data.scaled.l$value = as.numeric(as.character(data.scaled.l$value))

single.model.results = ddply(data.scaled.l, .var = c("variable"), .fun = function(x) {
  # K-fold validation for each model based on a single term 
  for(n in 1:folds.length) { 
    test.ids = folds$subsets[folds$which==n]
    test.trials = unique(data.scaled.l$trial)[test.ids]
    
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
  s.performance=mean(glmerperf)
  s.aic=mean(aic)
  s.coef=mean(coef)
  s.ref.sd=mean(ref.sd)
  print(x$variable[1])
  return(cbind(s.performance, s.aic, s.coef, s.ref.sd))
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
head(data.scaled.l)
dataPurged %>%
  ungroup() %>%
  tidyr::gather(key = variable, value = value, -trial, -success, -first_antigenicTypes_sd) -> data.scaled.l
data.scaled.l$value=as.numeric(data.scaled.l$value)


double.model.results = ddply(data.scaled.l, .var = c("variable"), .fun = function(x) {
  # K-fold validation for each model based on a single term 
  for(n in 1:folds.length) { 
    test.ids = folds$subsets[folds$which==n]
    test.trials = unique(data.scaled.l$trial)[test.ids]
    
    testdata = x[which(x$trial %in% test.trials), ]
    traindata = x[-which(x$trial %in% test.trials),]
    
    GLMER <- lme4::glmer(success ~ first_antigenicTypes_sd + value  + (1 | trial), data = traindata, family="binomial", 
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
  gather(key = variable, value = value, -trial, -success, -first_antigenicTypes_sd, -first_covBetaSigma_ratio) -> data.scaled.l
data.scaled.l$value=as.numeric(data.scaled.l$value)


third.model.results = ddply(data.scaled.l, .var = c("variable"), .fun = function(x) {
  # K-fold validation for each model based on a single term 
  for(n in 1:folds.length) { 
    test.ids = folds$subsets[folds$which==n]
    test.trials = unique(data.scaled.l$trial)[test.ids]
    
    testdata = x[which(x$trial %in% test.trials), ]
    traindata = x[-which(x$trial %in% test.trials),]
    
    GLMER <- lme4::glmer(success ~ first_antigenicTypes_sd + first_covBetaSigma_ratio + value  + (1 | trial), data = traindata, family="binomial", 
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
  gather(key = variable, value = value, -trial, -success, -selected.antigen, 
         -second_meanBeta_sd,
         -second_antigenicTypes_ratio,
         -first_diversity_sd,
         -second_covBetaSigma_ratio,
         -first_normalize.I_sd,
         -second_diversity_ratio,
         -first_meanBeta_ratio,
         -first_varBeta_diff) -> data.scaled.l
data.scaled.l$value=as.numeric(data.scaled.l$value)
# Remember to change the coefficien

nine.model.results = ddply(data.scaled.l, .var = c("variable"), .fun = function(x) {
  # K-fold validation for each model based on a single term 
  for(n in 1:folds.length) { 
    test.ids = folds$subsets[folds$which==n]
    test.trials = unique(data.scaled.l$trial)[test.ids]
    
    testdata = x[which(x$trial %in% test.trials), ]
    traindata = x[-which(x$trial %in% test.trials),]
    
     GLMER <- lme4::glmer(success ~  
                          second_meanBeta_sd +
                          second_antigenicTypes_ratio +
                          first_diversity_sd +
                          second_covBetaSigma_ratio +
                          first_normalize.I_sd+
                          second_diversity_ratio+
                          first_meanBeta_ratio+
                          first_varBeta_diff + 
                          value + (1 | trial), data = traindata, family="binomial", 
                           control = glmerControl(optimizer="bobyqa"),nAGQ=10)
    #GLMER <- glm(success ~ diversity + antigenicTypes + covBetaSigma + value, data = traindata, family="binomial")
    
    aic[n] <- extractAIC(GLMER)[2]
    coef[n] <- summary(GLMER)$coefficients[10,4]
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

nine.model.results

plot_pred_type_distribution <- function(df, threshold) {
  v <- rep(NA, nrow(df))
  v <- ifelse(df$pred >= threshold & df$survived == 1, "TP", v)
  v <- ifelse(df$pred >= threshold & df$survived == 0, "FP", v)
  v <- ifelse(df$pred < threshold & df$survived == 1, "FN", v)
  v <- ifelse(df$pred < threshold & df$survived == 0, "TN", v)
  
  df$pred_type <- v
  
  ggplot(data=df, aes(x=survived, y=pred)) + 
    geom_violin(fill=rgb(1,1,1,alpha=0.6), color=NA) + 
    geom_jitter(aes(color=pred_type), alpha=0.6) +
    geom_hline(yintercept=threshold, color="red", alpha=0.6) +
    scale_color_discrete(name = "type") +
    labs(title=sprintf("Threshold at %.2f", threshold))
}


#################### anova tests
dataPurged %>%
  dplyr::select(trial, success, selected.antigen, 
                second_meanBeta_sd, 
                second_antigenicTypes_ratio,
                first_diversity_sd,
                second_covBetaSigma_ratio,
                first_normalize.I_sd,
                second_diversity_ratio,
                first_meanBeta_ratio,
                first_varBeta_diff) -> data.scaled.l


for(n in 1:folds.length) { 
  test.ids = folds$subsets[folds$which==n]
  test.trials = unique(data.scaled.l$trial)[test.ids]
  
  #  testdata = x[which(x$trial %in% test.trials), ]
  #  traindata = x[-which(x$trial %in% test.trials),]
  
  testdata = data.scaled.l[which(data.scaled.l$trial %in% test.trials), ]
  traindata = data.scaled.l[-which(data.scaled.l$trial %in% test.trials),]
  traindata = na.omit(traindata)  

  model1 <- lme4::glmer(success ~  second_meanBeta_sd +
                                    second_antigenicTypes_ratio +
                                    first_diversity_sd +
                                    second_covBetaSigma_ratio +
                                    first_normalize.I_sd+
                                    second_diversity_ratio+
                                    first_meanBeta_ratio +
                                    first_varBeta_diff + (1 | trial), data = traindata, family="binomial", 
                                  control = glmerControl(optimizer="bobyqa"),nAGQ=10)
  
  model2 <-  lme4::glmer(success ~  second_meanBeta_sd +
                                                    second_antigenicTypes_ratio +
                                                    first_diversity_sd +
                                                    second_covBetaSigma_ratio +
                                                    first_normalize.I_sd+
                                                    second_diversity_ratio+
                                                    first_meanBeta_ratio+
                                                    first_varBeta_diff + 
                                                    first_meanBeta_ratio*first_varBeta_diff + (1 | trial), data = traindata, family="binomial", 
                                                  control = glmerControl(optimizer="bobyqa"),nAGQ=10)
  
  print(anova(model1, model2, test = "Chisq"))
}


summary(model2)

######################### Interaction Results 
dataPurged %>%
  select(trial, success, second_meanBeta_sd, 
         first_diversity_sd, first_antigenicTypes_ratio) -> data.scaled.l

coeff=rep(NA,folds.length)

for(n in 1:folds.length) { 
    test.ids = folds$subsets[folds$which==n]
    test.trials = unique(data.scaled.l$trial)[test.ids]
    
    testdata = data.scaled.l[which(data.scaled.l$trial %in% test.trials), ]
    traindata = data.scaled.l[-which(data.scaled.l$trial %in% test.trials),]

    model.null <- lme4::glmer(success ~ first_antigenicTypes_sd + first_covBetaSigma_ratio + first_diversity_sd +  (1 | trial),
                          data = traindata, family="binomial", 
                          control = glmerControl(optimizer="bobyqa"),nAGQ=10)

    model.int <- lme4::glmer(success ~ first_antigenicTypes_sd + first_covBetaSigma_ratio + first_diversity_sd + 
                              first_covBetaSigma_ratio*first_diversity_sd + (1 | trial),
                              data = traindata, family="binomial", 
                              control = glmerControl(optimizer="bobyqa"),nAGQ=10)
    print(summary(model.null))
  # coeff[n,] <- summary(model.9)[[10]][9:11,4]
    coef[n] <- summary(model.int)$coefficients[5,4]
    ref.sd[n] <- sqrt(VarCorr(model.int)$trial[1])
    aic[n] <- extractAIC(model.int)[2]
    glmer.probs <- predict(model.int, newdata=testdata, type="response", allow.new.levels=TRUE)
    glmer.ROC <- roc(predictor=glmer.probs, response=testdata$success)
    glmerperf[n] <- glmer.ROC$auc

   # print(anova(model.null, model.int, test="Chisq"))
  }

mean(glmerperf)
mean(aic)
mean(coef)
mean(ref.sd)
#colMeans(coeff)
               

plot(glmer.ROC, 0.7, 1,2)
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
y.offset = .4)http://127.0.0.1:25378/graphics/plot_zoom_png?width=1200&height=707
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


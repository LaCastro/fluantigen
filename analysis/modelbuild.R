library(plyr)
library(dplyr)
library(tidyverse)
library(broom)
library(lme4)
library(pROC)
library(caret)
library(cvTools)
library(cowplot)

#########################################################################################
predictorNames = colnames(freq.first.check)
exclude.emerge = c("freq", "day", "cases", "dominant.type", "cumulativeTypes", "antigentype","N",
                   "totalS", "totalI", "I", "R", "totalCases", "simDay")
                   
data = freq.first.check[, -which(predictorNames%in%exclude.emerge)]


data <- within(data, {
  success <- factor(success, levels=c("Est.", "Transient"), labels = c(1,0))
  name <- as.factor(name)
})

# Take our any data that's before 3 years
dataScaled = data
excluded.variables = c("netau")
dataScaled = dataScaled[,-which(colnames(dataScaled) %in% excluded.variables)]

factor.variables = which(sapply(dataScaled,is.factor)==TRUE)
#dataScaled$netau[!is.finite(dataScaled$netau)] <- NA
dataScaled[,-factor.variables] <- lapply(dataScaled[,-factor.variables],scale)
str(dataScaled)
trial.all <- glmer(success ~.-name + (1| name), data = dataScaled, family = binomial,
                   control = glmerControl(optimizer="bobyqa"),nAGQ=0)
vif.results = data.frame(vif.mer(trial.all))

vif.excluded = c("individual.meanMut", "S", "meanBeta",
                 "individual.varMut", "individual.meanSigma",
                 "individual.meanBeta","meanLoad",
                 "varSigma", "individual.meanR", "meanSigma",
                 "varR", "individual.varR",
                 "individual.varBeta", "meanR")

dataPurged = dataScaled[,-which(colnames(dataScaled) %in% vif.excluded)]
trial.purged <- lme4::glmer(success ~.-name + (1 | name), data = dataPurged, family = binomial,
                   control = glmerControl(optimizer="bobyqa"),nAGQ=0)
vif.results = data.frame(vif.mer(trial.purged))

################### Doing K-cross validation 
outcomeName = "success"

set.seed(62417)
folds.length = 5;
folds = cvFolds(n = length(unique(dataPurged$name)), K = folds.length, type = "random")

glmerperf=rep(0, folds.length); glmperf=glmerperf;
data.summary.mat = data.frame(matrix(nrow = folds.length, ncol  = 10))
colnames(data.summary.mat) = colnames(data.summary)

dataPurged %>%
  ungroup() %>%
  gather(key = variable, value = value, -success, -name) -> data.scaled.l
data.scaled.l$value = as.numeric(as.character(data.scaled.l$value))

single.model.results = ddply(data.scaled.l, .var = c("variable"), .fun = function(x) {
  # K-fold validation for each model based on a single term 
  for(n in 1:folds.length) { 
    test.ids = folds$subsets[folds$which==n]
    test.trials = unique(data.scaled.l$name)[test.ids]
    testdata = x[which(x$name %in% test.trials), ]
    traindata = x[-which(x$name %in% test.trials),]
    
    GLMER <- lme4::glmer(success ~ value  + (1 | name), data = traindata, family="binomial", 
                         control = glmerControl(optimizer="bobyqa"),nAGQ=0)
   
    data.summary1 = glance(GLMER); data.summary2 = tidy(GLMER)[2,2:5]
    data.summary = bind_cols(data.summary1,data.summary2)
    data.summary.mat[n,] = data.summary
    glmer.probs <- predict(GLMER, newdata=testdata, type="response", allow.new.levels=TRUE)
    glmer.ROC <- roc(predictor=glmer.probs, response=testdata$success)
    glmerperf[n] <- glmer.ROC$auc
  }
  s.performance=mean(glmerperf)
  s.data.summary = colMeans(as.data.frame(data.summary.mat))
  print(x$variable[1])
  return.data = t(c(per = s.performance, s.data.summary))
  return(return.data)
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
  
  
# model.emerge <- glm(success~distance+mutLoad + normalize.I + meanLoad+ meanLoad*mutLoad + mutLoad*normalize.I,
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
  ungroup() %>%
  tidyr::gather(key = variable, value = value,  -success, -name, -ratio.meanR) -> data.scaled.l
data.scaled.l$value=as.numeric(data.scaled.l$value)

double.model.results = ddply(data.scaled.l, .var = c("variable"), .fun = function(x) {
  # K-fold validation for each model based on a single term 
  for(n in 1:folds.length) { 
    test.ids = folds$subsets[folds$which==n]
    test.trials = unique(data.scaled.l$name)[test.ids]
    testdata = x[which(x$name %in% test.trials), ]
    traindata = x[-which(x$name %in% test.trials),]
    
    GLMER <- lme4::glmer(success ~ value  + ratio.meanR + (1 | name),
                         data = traindata, family="binomial", 
                         control = glmerControl(optimizer="bobyqa"),nAGQ=0)
    
    data.summary1 = glance(GLMER); data.summary2 = tidy(GLMER)[2,2:5]
    data.summary = bind_cols(data.summary1,data.summary2)
    data.summary.mat[n,] = data.summary
    glmer.probs <- predict(GLMER, newdata=testdata, type="response", allow.new.levels=TRUE)
    glmer.ROC <- roc(predictor=glmer.probs, response=testdata$success)
    glmerperf[n] <- glmer.ROC$auc
  }
  performance=mean(glmerperf)
  data.summary = colMeans(as.data.frame(data.summary.mat))
  print(x$variable[1])
  return.data = t(c(per = performance, data.summary))
  return(return.data)
})
double.model.results %>% 
  arrange(desc(per))


############################### third term results
dataPurged %>%
  ungroup() %>%
  tidyr::gather(key = variable, value = value,  -success, -name, -ratio.meanR, -ratio.meanSigma) -> data.scaled.l
data.scaled.l$value=as.numeric(data.scaled.l$value)

third.model.results = ddply(data.scaled.l, .var = c("variable"), .fun = function(x) {
  # K-fold validation for each model based on a single term 
  for(n in 1:folds.length) { 
    test.ids = folds$subsets[folds$which==n]
    test.trials = unique(data.scaled.l$name)[test.ids]
    testdata = x[which(x$name %in% test.trials), ]
    traindata = x[-which(x$name %in% test.trials),]
    
    GLMER <- lme4::glmer(success ~ value  + ratio.meanSigma+ratio.meanR + (1 | name),
                         data = traindata, family="binomial", 
                         control = glmerControl(optimizer="bobyqa"),nAGQ=0)
    
    data.summary1 = glance(GLMER); data.summary2 = tidy(GLMER)[2,2:5]
    data.summary = bind_cols(data.summary1,data.summary2)
    data.summary.mat[n,] = data.summary
    glmer.probs <- predict(GLMER, newdata=testdata, type="response", allow.new.levels=TRUE)
    glmer.ROC <- roc(predictor=glmer.probs, response=testdata$success)
    glmerperf[n] <- glmer.ROC$auc
  }
  performance=mean(glmerperf)
  data.summary = colMeans(as.data.frame(data.summary.mat))
  print(x$variable[1])
  return.data = t(c(per = performance, data.summary))
  return(return.data)
})
third.model.results %>% 
  arrange(desc(per))

third.model.results

############################### Four and greater term results
dataPurged %>%
  ungroup() %>%
  tidyr::gather(key = variable, value = value,  -success, -name, 
                -ratio.meanR, -ratio.meanSigma, -ratio.mutation,
                -individual.varSigma, -infected, -ratio.varR) -> data.scaled.l
data.scaled.l$value=as.numeric(data.scaled.l$value)

seven.model.results = ddply(data.scaled.l, .var = c("variable"), .fun = function(x) {
  # K-fold validation for each model based on a single term 
  for(n in 1:folds.length) { 
    test.ids = folds$subsets[folds$which==n]
    test.trials = unique(data.scaled.l$name)[test.ids]
    testdata = x[which(x$name %in% test.trials), ]
    traindata = x[-which(x$name %in% test.trials),]
    
    GLMER <- lme4::glmer(success ~ value  + individual.varSigma + ratio.mutation + 
                           ratio.meanSigma+ratio.meanR + infected + ratio.varR +
                           ratio.varBeta + (1 | name),
                         data = traindata, family="binomial", 
                         control = glmerControl(optimizer="bobyqa"),nAGQ=0)
    
    data.summary1 = glance(GLMER); data.summary2 = tidy(GLMER)[2,2:5]
    data.summary = bind_cols(data.summary1,data.summary2)
    data.summary.mat[n,] = data.summary
    glmer.probs <- predict(GLMER, newdata=testdata, type="response", allow.new.levels=TRUE)
    glmer.ROC <- roc(predictor=glmer.probs, response=testdata$success)
    glmerperf[n] <- glmer.ROC$auc
  }
  performance=mean(glmerperf)
  data.summary = colMeans(as.data.frame(data.summary.mat))
  print(x$variable[1])
  return.data = t(c(per = performance, data.summary))
  return(return.data)
})
seven.model.results %>% 
  arrange(desc(per))


############################## Interaction
dataPurged %>%
  ungroup() %>%
  select(success, name, ratio.meanR, ratio.meanSigma, ratio.mutation,
                individual.varSigma, infected, ratio.varR) -> data.scaled.l

for(n in 1:folds.length) { 
  test.ids = folds$subsets[folds$which==n]
  test.trials = unique(data.scaled.l$name)[test.ids]
  
  testdata = data.scaled.l[which(data.scaled.l$name %in% test.trials), ]
  traindata = data.scaled.l[-which(data.scaled.l$name %in% test.trials),]
  
  model.null <- lme4::glmer(success ~ ratio.meanR + ratio.meanSigma + ratio.mutation +
                              individual.varSigma +infected + ratio.varR  +  (1 | name),
                            data = traindata, family="binomial", 
                            control = glmerControl(optimizer="bobyqa"),nAGQ=10)
  
  model.int <- lme4::glmer(success ~ ratio.meanR*ratio.meanSigma + 
                             ratio.meanR + ratio.meanSigma + ratio.mutation +
                             individual.varSigma +infected + ratio.varR  +  (1 | name),
                           data = traindata, family="binomial", 
                           control = glmerControl(optimizer="bobyqa"),nAGQ=10)
 
  anova.test = anova(model.null, model.int)
  data.summary1 = glance(model.int); data.summary2 = tidy(model.int)[8,2:5]
  data.summary = bind_cols(data.summary1,data.summary2)
  data.summary.mat[n,] = data.summary
  glmer.probs <- predict(model.int, newdata=testdata, type="response", allow.new.levels=TRUE)
  glmer.ROC <- roc(predictor=glmer.probs, response=testdata$success)
  glmerperf[n] <- glmer.ROC$auc
} 
performance=mean(glmerperf)

data.summary = colMeans(as.data.frame(data.summary.mat))
return.data = t(c(per = performance, data.summary, anova.test = anova.test$`Pr(>Chisq)`))
  return(return.data)
  
  return.data = t(c(per = performance, data.summary))
  
  # print(anova(model.null, model.int, test="Chisq"))
}









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


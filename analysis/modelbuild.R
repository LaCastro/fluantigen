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
freq.df.subset %>%
  filter(freq == "freq.01") -> freq.1.subset

predictorNames = colnames(freq.1.subset)

data = freq.1.subset

data <- within(data, {
  success <- factor(success ) #, levels=c("Est.", "Transient"))
  name <- as.factor(name)
  freq <- as.factor(freq)
  antigentype <- as.factor(antigentype)
})


dataScaled = data %>% select(-freq, -antigentype, -simDay, -day)
dataScaled$netau[is.infinite(dataScaled$netau)] <- NA

factor.variables = which(sapply(dataScaled,is.factor)==TRUE)
dataScaled[,-factor.variables] <- lapply(dataScaled[,-factor.variables],scale)


########## Fit first model 
trial.all <- glmer(success ~.-name + (1| name), data = dataScaled, family = binomial,
                   control = glmerControl(optimizer="bobyqa"),nAGQ=0)


vif.results = data.frame(vif.mer(trial.all))

vif.excluded = c("individual.meanMut", "individual.meanBeta", "individual.varR",
                 "individual.varMut", "totalS", "meanSigma", "varSigma" , "individual.meanR",
                 "meanR", "ratio.meanR", "individual.varBeta", "meanLoad", "ratio.varBeta",
                 "ratio.meanBeta","antigenicTypes", "totalI")



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
colnames(data.summary.mat) = c("sigma", "logLik", "AIC", "BIC", "deviance", "df.residual", "estimate", "std.error", "statistic", "p.value")

dataPurged$success <- relevel(dataPurged$success, ref = "Transient")
dataPurged %>%
  ungroup() %>%
  gather(key = variable, value = value, -success, -name, -netau) -> data.scaled.l

data.scaled.l$value = as.numeric(as.character(data.scaled.l$value))

full.variables = unique(data.scaled.l$variable)

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
variable.1 = single.model.results %>% arrange(AIC) %>% slice(1)

variable.list = c(variable.1$variable)
variable.per = c(variable.1$per)
variable.aic = c(variable.1$AIC)

########################### second.model.results 


variable.list

dataPurged %>%
  ungroup() %>%
  tidyr::gather(key = variable, value = value,  -success, -name,
                -netau, -ratio.meanSigma,-ratio.mutation) -> data.scaled.l
data.scaled.l$value=as.numeric(data.scaled.l$value)

variable.list = c(variable.list, "ratio.mutation")

fixed.part.1 = "success ~ value + "
variable.part = paste(variable.list, collapse = "+")
fixed.part.2 = " + (1 | name)"

formula = paste0(fixed.part.1, variable.part, fixed.part.2)


double.model.results = ddply(data.scaled.l, .var = c("variable"), .fun = function(x) {
  # K-fold validation for each model based on a single term 
  for(n in 1:folds.length) { 
    test.ids = folds$subsets[folds$which==n]
    test.trials = unique(data.scaled.l$name)[test.ids]
    testdata = x[which(x$name %in% test.trials), ]
    traindata = x[-which(x$name %in% test.trials),]
    
    GLMER <- lme4::glmer(formula,
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
  arrange(AIC)


############################### third term results
dataPurged %>%
  ungroup() %>%
  tidyr::gather(key = variable, value = value,-success, -name, 
                -netau, -ratio.meanSigma, -ratio.mutation) -> data.scaled.l
data.scaled.l$value=as.numeric(data.scaled.l$value)

third.model.results = ddply(data.scaled.l, .var = c("variable"), .fun = function(x) {
  # K-fold validation for each model based on a single term 
  for(n in 1:folds.length) { 
    test.ids = folds$subsets[folds$which==n]
    test.trials = unique(data.scaled.l$name)[test.ids]
    testdata = x[which(x$name %in% test.trials), ]
    traindata = x[-which(x$name %in% test.trials),]
    
    GLMER <- lme4::glmer(success ~ value  + 
                           ratio.mutation + ratio.meanSigma + (1 | name),
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
  arrange(AIC)

############################### Four and greater term results
dataPurged %>%
  ungroup() %>%
  tidyr::gather(key = variable, value = value,  -success, -name, 
              -ratio.meanSigma, -ratio.mutation, -varR, -netau, -individual.meanSigma,
              -individual.varSigma, -ratio.varR, -meanBeta, -infected) -> data.scaled.l

data.scaled.l$value=as.numeric(data.scaled.l$value)
unique(data.scaled.l$variable)

nine.results = ddply(data.scaled.l, .var = c("variable"), .fun = function(x) {
  # K-fold validation for each model based on a single term 
  for(n in 1:folds.length) { 
    test.ids = folds$subsets[folds$which==n]
    test.trials = unique(data.scaled.l$name)[test.ids]
    testdata = x[which(x$name %in% test.trials), ]
    traindata = x[-which(x$name %in% test.trials),]
    
    GLMER <- lme4::glmer(success ~ value  + 
                           ratio.meanSigma + 
                           ratio.mutation + varR +
                           individual.meanSigma + individual.varSigma + 
                           ratio.varR + meanBeta + infected + (1 | name),
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
  return(t(c(per = performance, data.summary)))
})
nine.model.results %>% 
  arrange(AIC)



############################## Anova 
dataPurged %>%
  ungroup() %>%
  dplyr::select(success, name,
                netau, ratio.meanSigma, ratio.mutation, varR,
                individual.meanSigma, individual.varSigma, ratio.varR,
                infected, entropy, diversity) -> data.scaled.l
anova.results = rep(0, folds.length)

for(n in 1:folds.length) { 
    test.ids = folds$subsets[folds$which==n]
    test.trials = unique(data.scaled.l$name)[test.ids]
    testdata = data.scaled.l[which(data.scaled.l$name %in% test.trials), ]
    traindata = data.scaled.l[-which(data.scaled.l$name %in% test.trials),]
    
    GLMER.null <- lme4::glmer(success ~ ratio.meanSigma + ratio.mutation + varR + 
                                individual.meanSigma + individual.varSigma + 
                                ratio.varR + infected + entropy + (1 | name),
                         data = traindata, family="binomial", 
                         control = glmerControl(optimizer="bobyqa"),nAGQ=0)
    GLMER.test <- lme4::glmer(success ~ ratio.meanSigma + ratio.mutation + varR + 
                                individual.meanSigma + individual.varSigma + 
                                ratio.varR + infected + entropy + diversity +  (1 | name),
                              data = traindata, family="binomial", 
                              control = glmerControl(optimizer="bobyqa"),nAGQ=0)
  anova.test = anova(GLMER.null, GLMER.test)
  anova.results[n] = anova.test$`Pr(>Chisq)`[2]
  }
mean(anova.results)

############################## Interaction
dataPurged %>%
  ungroup() %>%
  dplyr::select(success, name,
                ratio.meanSigma, ratio.mutation, varR, individual.meanSigma, 
                individual.varSigma, ratio.varR, meanBeta, infected) -> data.scaled.l
anova.test=rep(0,folds.length)
for(n in 1:folds.length) { 
  test.ids = folds$subsets[folds$which==n]
  test.trials = unique(data.scaled.l$name)[test.ids]
  testdata = data.scaled.l[which(data.scaled.l$name %in% test.trials), ]
  traindata = data.scaled.l[-which(data.scaled.l$name %in% test.trials),]
  
  model.null <- lme4::glmer(success ~ individual.varSigma*ratio.varR + 
                              varR*individual.meanSigma + 
                              ratio.mutation*ratio.varR +
                              ratio.mutation*meanBeta + 
                              ratio.meanSigma + ratio.mutation + varR +
                              individual.meanSigma + individual.varSigma + 
                              ratio.varR + meanBeta + infected + (1 | name),
                            data = traindata, family="binomial", 
                            control = glmerControl(optimizer="bobyqa"),nAGQ=0)
  model.int <- lme4::glmer(success ~meanBeta*infected+ 
                             ratio.mutation * meanBeta + ratio.mutation*ratio.varR + 
                             varR*individual.meanSigma + individual.varSigma*ratio.varR  +
                             ratio.meanSigma + ratio.mutation + varR + individual.meanSigma + 
                             individual.varSigma +  ratio.varR + meanBeta + infected + (1 | name),
                           data = traindata, family="binomial", 
                           control = glmerControl(optimizer="bobyqa"),nAGQ=0)
  model.test = anova(model.null, model.int)
  anova.test[n] = model.test$`Pr(>Chisq)`[2]
  data.summary1 = glance(model.int); data.summary2 = tidy(model.int)[10,2:5]
  data.summary.mat[n,] = bind_cols(data.summary1,data.summary2)
  glmer.probs <- predict(model.int, newdata=testdata, type="response", allow.new.levels=TRUE)
  glmer.ROC <- roc(predictor=glmer.probs, response=testdata$success)
  glmerperf[n] <- glmer.ROC$auc
} 
mean(anova.test)
mean(glmerperf)
colMeans(as.data.frame(data.summary.mat))


summary(model.int)





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

freq.explore %>%
  filter(freq == ".01") -> freq.01

levels(dataScaled$success)

dataScaled$success <- relevel(dataScaled$success, ref = "Transient")
freq.01.model <- lme4::glmer(success~ratio.meanR + ratio.meanSigma+ratio.mutation+individual.varSigma+infected+ratio.varR +
                               ratio.varBeta + infected*ratio.varBeta+ratio.meanR*ratio.mutation + infected*ratio.varR + (1 | name),
                             data = dataScaled, family = "binomial", control = glmerControl(optimizer="bobyqa", nAGQ=10))

dataScaled$success <- relevel(dataScaled$success, ref = "Transient")
freq.02.model <- lme4::glmer(success~ ratio.meanR + ratio.meanSigma + ratio.mutation + 
                               individual.varSigma + antigenicTypes+ratio.varR + 
                               tmrca + ratio.meanR*ratio.mutation + antigenicTypes*ratio.varR + (1 | name),
                             data = dataScaled, family = "binomial", control = glmerControl(optimizer="bobyqa", nAGQ=10))
dataScaled$success <- relevel(dataScaled$success, ref = "Transient")
freq.05.model <- lme4::glmer(success~ratio.meanR + ratio.meanSigma+ratio.mutation+ratio.varSigma +
                               ratio.varR + antigenicTypes +  (1 | name), data = dataScaled, family = "binomial", 
                                                         control = glmerControl(optimizer="bobyqa", nAGQ=10))


dataScaled$success <- relevel(dataScaled$success, ref = "Transient")
colnames(dataScaled)
diff.second.model <- lme4::glmer(success~ratio.meanBeta+covBetaSigma+dominant.freq+varR+individual.varSigma+infected + (1 | name),
                                 data = dataScaled, family = "binomial", control = glmerControl(optimizer="bobyqa", nAGQ=10))


dataScaled$success <- relevel(dataScaled$success, ref = "Transient")
three.point.model <- lme4::glmer(success~five_ratio.meanR + five_meanR + five_varR + second_meanR + 
                                  five_individual.varSigma+five_ratio.varBeta + first_ratio.meanSigma+five_ratio.mutation +
                                  five_ratio.meanR*second_meanR + (1 | name), data = dataScaled, family = "binomial",
                                 control = glmerControl(optimizer="bobyqa", nAGQ=10)) 

second.point.model <- lme4::glmer(success~first_ratio.meanSigma + first_ratio.mutation + second_infected + 
                                    first_varR + first_meanR + first_varR*first_meanR + (1 | name), 
                                  data = dataScaled, family = "binomial", control = glmerControl(optimizer="bobyqa", nAGQ=10))

dataScaled$success <- relevel(dataScaled$success, ref = "Transient")
diff.two.point.model <- lme4::glmer(success~second.diff.ratio.meanBeta+ first.diff.varR + second.diff.dominant.freq + first.diff.ratio.meanBeta + 
                                      second.diff.I + second.diff.individual.varSigma+second.diff.covBetaSigma+second.diff.dominant.freq*second.diff.I + (1|name),
                                    data = dataScaled, family = "binomial", control = glmerControl(optimizer="bobyqa", nAGQ=10))

dataScaled$success <- relevel(dataScaled$success, ref = "Transient")
diff.second.growth.model <- lme4::glmer(success~ratio.meanBeta + covBetaSigma+dominant.freq+varR +individual.varSigma +
                                          infected + (1 | name), data = dataScaled, family = "binomial", control = glmerControl(optimizer="bobyqa", nAGQ=10))



dataScaled$success <- relevel(dataScaled$success, ref = "Transient")
diff.first.growth.model <- lme4::glmer(success~ varR + ratio.meanBeta + individual.meanSigma + covBetaSigma + (1 | name), data = dataScaled, family = "binomial", control = glmerControl(optimizer="bobyqa", nAGQ=10))



summary(freq.05.model)
summary(freq.02.model)
summary(three.point.model)
summary(second.point.model)
summary(diff.two.point.model)
summary(diff.second.growth.model)
summary(diff.first.growth.model)

library(sjPlot)
library(sjmisc)
library(sjlabelled)

sjp.glmer(freq.05.model, type = "fe") 
sjp.glmer(freq.02.model, type = "fe")
sjp.glmer(three.point.model, type = "fe")
sjp.glmer(second.point.model, type = "fe")
sjp.glmer(diff.two.point.model, type = "fe")  
sjp.glmer(diff.first.growth.model, type = "fe")
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




##########################################

model.summary <- read.csv("../data/model.csv", header = TRUE)

model.types = c("difference.first", "difference.second", "diff.first+second",
                "state.01", "state.02", "state.03", "state.01+02+05", "state.01+02")
model.summary %>%
  ggplot(aes(x=as.factor(Term), y = Per, group = Model, color = Model)) + geom_line(size = 1.2) +
  labs(x = "Number of Terms", y = "5-Fold CV AUC") + scale_color_brewer(palette = "Spectral", labels = model.types) -> model.summary.plot

save_plot(fig = model.summary.plot, filename = "exploratory.figs/model.summary.plot.pdf", base_aspect_ratio = 1.6)

freq.models = c("freq.01", "freq.02", "freq.05")
diff.models = c("diff.first", "diff.second")
time.point.freq = c("two.point", "three.point")

unique(model.summary$Model)

model.summary %>%
  filter(Model== "three.point") %>%
  arrange(Term)



  
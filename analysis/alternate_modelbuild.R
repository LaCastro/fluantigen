library(plyr)
library(dplyr)
library(tidyverse)
library(lme4)
library(pROC)
library(caret)
library(cvTools)
library(cowplot)

#########################################################################################
predictorNames = colnames(meta.data)
exclude.emerge = c("oriAntigen", "postAntigen", "cases", "cumulativeTypes", 
                   "final.max", "dominant.type", "life.length", "max.I", "min.I",
                   "days.above", "infected", "N", "R")

surv.list = list(two = freq.two, five = freq.five)

data = antigen.data[, -which(predictorNames%in%exclude.emerge)]

excluded = c("N","R","I", "totalS", "totalI", "totalCases", "cases", "totalCases", "cumulativeTypes", "dominant.type") 
freq = freq.two[, -which(colnames(freq.two)%in%excluded)]


freq <- within(freq, {
  antigentype <- factor(antigentype)
  success <- factor(success, levels=c("Est.", "Transient"), labels = c(1,0))
  name <- factor(name)
  day <- as.numeric(as.character(day))
  #  mutLoad <- as.numeric(as.character(mutLoad))
  #  distance <- as.numeric(as.character(distance))
  #  simDay <- as.numeric(as.character(simDay))
})

freq %>%
  filter(day > (3*365)) -> freq


##### When the simDay isn't in terms of the burn in day, don't need the -1 and -3 
#dataBurn %>%
#  mutate(time.of.year = (simDay)-floor((x = day/365))) %>%
#  mutate(quarter = ifelse(time.of.year < .25, 1,
#                          ifelse(time.of.year < .5, 2,
#                                 ifelse(time.of.year < .75, 3,4)))) -> dataScaled

#dataScaled$quarter = as.factor(dataScaled$quarter)

# Take our any data that's before 3 years
dataScaled = freq

#dataScaled%>%
#  group_by(name)%>%
#  summarize(mutations = n())%>%
#  summarize(mean.mutations = mean(mutations)/number.years.rep)

## Variables to remove -- netau and mean Load ratio 
factor.variables = which(sapply(dataScaled,is.factor)==TRUE)
dataScaled$netau[!is.finite(dataScaled$netau)] <- NA
dataScaled[,-factor.variables] <- lapply(dataScaled[,-factor.variables],scale)

trial.all <- glmer(success ~.-name-antigentype-day  + (1|name), data = dataScaled, family = binomial,
                   control = glmerControl(optimizer="bobyqa"),nAGQ=0)
vif.results = data.frame(vif.mer(trial.all)) 


vif.excluded = c("S", "meanBeta", "varBeta", "meanSigma")
dataPurged = dataScaled[,-which(colnames(dataScaled) %in% vif.excluded)]

trial.purged <- lme4::glmer(success ~.-name-antigentype-day  + (1 | name), data = dataPurged, family = binomial,
                            control = glmerControl(optimizer="bobyqa"),nAGQ=0)
data.frame(vif.mer(trial.purged))


################### Doing K-cross validation 
outcomeName = "success"

set.seed(62417)
folds.length = 5;
folds = cvFolds(n = length(unique(dataPurged$name)), K = folds.length, type = "random")

## Create Folds DV
folds.df.test = data.frame(cbind(fold = folds$which, folds$subsets)); colnames(folds.df.test)[2] = "row.index"


find_first_term <- function(freq, folds) { 

  ## Create Model List
  variables = colnames(dataPurged); variables = variables[-(which(variables %in% c("success", "name", "antigentype")))]
  model.formulas = as.list(paste0("success", "~", variables, "+ (1|name)"))

  
# Calculating Single Term   
single.term = map(model.formulas, five_fold_cv) %>%
  bind_rows() %>%
  arrange(p.value)
  
# Select First Model Term 
first.term = single.term %>% slice(1)

model.terms = first.term$term
### Generate Model List

variables = variables[-(which(variables %in% model.terms))]

model.formulas = as.list(paste0("success", "~", ))

generate_sig_terms_formula <- function(model.terms) {
  string = ""
  for(i in 1:length(model.terms) {
    string.term = model.terms[i]
    string = cbind(string, string.term)
    if(length(model.terms) > (i +1) {
      string = c(string, "+")
    })
  }
}




## Need function to fit based on folds 
five_fold_cv <- function(formula, folds.df = folds.df.test, data = dataPurged) {  

 fold.results = ddply(folds.df, .variables = "fold", function(fold) fit_single_fold(formula = formula, fold.n = fold, data = dataPurged)) 
 fold.results %>%
   summarize(term = term[1],
             p.value = mean(p.value),
             AIC = mean(AIC),
             auc = mean(auc),
             estimate = mean(estimate),
             std.error = mean(std.error)) 
}
fit_single_fold <- function(formula, fold.n, data = dataPurged) {
  testdata =  set_test_data(data, fold.n)
  traindata = set_train_data(data,fold.n)  
  
  GLMER <- lme4::glmer(formula, data = traindata, family="binomial", 
                       control = glmerControl(optimizer="bobyqa"),nAGQ=0)
  glmer.probs <- predict(GLMER, newdata=testdata, type="response", allow.new.levels=TRUE)
  glmer.ROC <- roc(predictor=glmer.probs, response=testdata$success)
  
  glance.data = glance(GLMER)
  tidy.data = tidy(GLMER)
  
  tidy.data %>%
    filter(!(term %in% c("(Intercept)", "sd_(Intercept).name"))) %>%
    slice(n()) %>% ### this will be the last one in the data frame
    bind_cols(glance.data) %>%
    bind_cols(auc = glmer.ROC$auc)
}
set_test_data <- function(data, fold) {
  test.ids = fold$row.index
  testnames = unique(data$name)[test.ids]
  testdata = data[which(data$name %in% testnames),]
  return(testdata)
}
set_train_data <- function(data, fold) {
  test.ids = fold$row.index
  testnames = unique(data$name)[test.ids]
  traindata = data[-which(data$name %in% testnames),]
  return(traindata)
}


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



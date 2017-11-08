library(plyr)
library(dplyr)
library(tidyverse)
library(broom)
library(lme4)
library(pROC)
library(caret)
library(cvTools)
library(cowplot)
#################################################################
population.level = c("totalS","infected", "antigenicTypes", "entropy", 
                     "diversity","antigenicDiversity", "tmrca", "dominant.freq",
                     "serialInterval")
antigen.level = c("individual.meanMut", "individual.meanBeta",
                  "individual.varBeta", "individual.varMut","individual.meanR", 
                  "individual.meanSigma", "individual.varR", "individual.varSigma")  
viral.popfitness = c("meanR", "varR", "meanSigma", "varSigma", "covBetaSigma", 
                     "meanBeta", "varBeta", "meanLoad")
ratio.fitness = c("ratio.mutation", "ratio.meanR", "ratio.varR", "ratio.meanBeta",
                  "ratio.varBeta", "ratio.meanSigma", "ratio.varSigma")


#########################################################################################
freq.df.subset %>%
  filter(freq == "freq.01") -> freq.1.subset

data = accelerate

data  <- within(data, {
  success <- factor(success) #, levels=c("Est.", "Transient"))
  name <- as.factor(name)
#  freq <- as.factor(freq)
  antigentype <- as.factor(antigentype)
})

#dataScaled = data %>% select(success, name, antigentype) # Snap shot 
dataScaled = data %>% select(-simDay,-totalI, -netau) # Growth Difference 

#dataScaled$netau[is.infinite(dataScaled$netau)] <- NA
factor.variables = which(sapply(dataScaled,is.factor)==TRUE)
dataScaled[,-factor.variables] <- lapply(dataScaled[,-factor.variables],scale)

########## Fit first model 
trial.all <- glmer(success ~.-name-antigentype + (1| name), data = dataScaled, family = binomial,
                   control = glmerControl(optimizer="bobyqa"),nAGQ=0)
vif.results = data.frame(vif.mer(trial.all))
# Snap Shot 
#vif.excluded = c("individual.meanMut", "individual.meanBeta", "individual.varR",
#                 "individual.varMut", "totalS", "meanSigma", "varSigma" , "individual.meanR",
#                 "meanR", "ratio.meanSigma", "individual.varBeta", "meanLoad", "ratio.varBeta",
#                 "ratio.meanBeta","antigenicTypes", "totalI", "day")
# Growth Differences 
vif.excluded = c("individual.meanMut", "individual.meanBeta", "individual.varR",
                 "individual.varMut", "totalS", "meanSigma", "varSigma" , "individual.meanR",
                 "meanR", "ratio.meanSigma", "individual.varBeta", "meanLoad", "ratio.varBeta",
                 "ratio.meanBeta","totalI")

#dataPurged = dataScaled

dataPurged = dataScaled[,-which(colnames(dataScaled) %in% vif.excluded)]
trial.purged <- lme4::glmer(success ~.-name-antigentype + (1 | name), data = dataPurged, family = binomial,
                   control = glmerControl(optimizer="bobyqa"),nAGQ=0)
vif.results = data.frame(vif.mer(trial.purged))
vif.results %>% arrange(desc(vif.mer.trial.purged.)) %>% head()

colnames(dataScaled)[which(!colnames(dataScaled) %in% colnames(dataPurged))]

################# Combing multiple time point data sets 
factor.variables = which(sapply(dataPurged,is.factor)==TRUE)
colnames(dataPurged)[-factor.variables] = paste0("accel.", colnames(dataPurged)[-factor.variables])

freq.1 = dataPurged
freq.3 = dataPurged
freq.5 = dataPurged
diff.1 = dataPurged
diff.2 = dataPurged
accel = dataPurged

freq.1 %>% gather(key = variable, value = value, -success,-name, -antigentype) -> freq.1.l
freq.3 %>% gather(key = variable, value = value, -success,-name, -antigentype) -> freq.3.l
freq.5 %>% gather(key = variable, value = value, -success,-name, -antigentype) -> freq.5.l
diff.2 %>% gather(key = variable, value = value, -success, -name, -antigentype) -> dataPurged.2
diff.1 %>% gather(key = variable, value = value, -success, -name, -antigentype)  -> dataPurged.1
accel %>% gather(key = variable, value = value, -success, -name, -antigentype) -> accel.l

the.whole.data = bind_rows(freq.1.l, freq.3.l, freq.5.l, dataPurged.1, dataPurged.2, accel.l) %>%
  spread(key = variable, value = value) %>% 
  select(-antigentype)

trial.multiple <- lme4::glmer(success ~.-name + (1 | name), data = the.whole.data, family = binomial,
                            control = glmerControl(optimizer="bobyqa"),nAGQ=0)
vif.results = data.frame(vif.mer(trial.multiple))

vif.excluded = c("freq.1individual.varSigma", "freq.5individual.varSigma", "freq.3meanBeta", "freq.5meanBeta",
                 "freq.3antigenicDiversity", "freq.3diversity", "freq.3individual.meanSigma", "freq.3entropy", 
                 "freq.3dominant.freq", "freq.3infected", "freq.3varBeta", "freq.3covBetaSigma", "freq.3ratio.mutation",
                 "freq.5antigenicDiversity", "freq.5diversity", "freq.1individual.meanSigma", "gp.1antigenicTypes", 
                 "freq.3ratio.meanR", "gp.1individual.meanSigma", "freq.1entropy", "freq.3tmrca", "freq.5individual.meanSigma",
                 "freq.1covBetaSigma", "freq.5covBetaSigma", "freq.5dominant.freq", "gp.1dominant.freq", "freq.3varR")
the.whole.data.purged = the.whole.data[,-which(colnames(the.whole.data) %in% vif.excluded)]
trial.purged <- lme4::glmer(success ~.-name + (1 | name), data = the.whole.data.purged, family = binomial,
                            control = glmerControl(optimizer="bobyqa"),nAGQ=0)
vif.results = data.frame(vif.mer(trial.purged))


dataPurged = the.whole.data.purged

included.variables = rownames(vif.results)
write.csv(included.variables, file = "../results/included.variables.csv")






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






  
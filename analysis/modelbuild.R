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
full.level = c("totalS", "infected", "antigenicTypes", "entropy", "diversity", "antigenicDiversity", "tmrca",
               "dominant.freq", "serialInterval", "meanR")
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

#########################################################################
freq.df.subset %>%
  filter(freq == "freq.02") -> freq.2.subset

data = accelerate

data  <- within(data, {
  success <- factor(success) #, levels=c("Est.", "Transient"))
  name <- as.factor(name)
#  freq <- as.factor(freq)
  antigentype <- as.factor(antigentype)
})
 
dataScaled = data %>% select(-simDay, -netau, -infected) %>%
  mutate_if(is.numeric, scale)# Snap shot 
dataPurged = dataScaled

########## Fit first model 
#trial.all <- glmer(success ~.-name-antigentype + (1| name), data = dataScaled, family = binomial,
#                   control = glmerControl(optimizer="bobyqa"),nAGQ=0)
#vif.results = data.frame(vif.mer(trial.all))
#View(vif.results)

# Snap Shot 
#vif.excluded = c("individual.meanMut", "individual.varMut", "individual.meanBeta", "individual.meanR", "individual.meanSigma",
#                 "ratio.meanSigma", "meanSigma", "covBetaSigma", "meanBeta", "individual.varR", "individual.varBeta", "antigenicTypes",
#                 "ratio.varR")

# 1/22/
#vif.excluded = c("individual.meanMut", "individual.varMut", "individual.meanBeta", "individual.meanR", "individual.meanSigma",
#                 "ratio.meanSigma", "meanSigma", "covBetaSigma", "meanBeta", "individual.varR", "individual.varBeta", "ratio.varR")

# 1/22/ Growth 
#vif.excluded = c("individual.meanMut", "individual.varBeta", "individual.varR", "ratio.meanR", "meanR", "ratio.meanBeta",
#                 "individual.meanR", "ratio.meanSigma", "varSigma", "individual.varMut")

#dataPurged = dataScaled[,-which(colnames(dataScaled) %in% vif.excluded)]
#trial.purged <- lme4::glmer(success ~.-name-antigentype + (1 | name), data = dataPurged, family = binomial,
#                   control = glmerControl(optimizer="bobyqa"),nAGQ=0)
#vif.results = data.frame(vif.mer(trial.purged))
#View(vif.results)
#vif.results %>% arrange(desc(vif.mer.trial.purged.)) %>% head()
#colnames(dataScaled)[which(!colnames(dataScaled) %in% colnames(dataPurged))]

################# Combing multiple time point data sets 
factor.variables = which(sapply(dataPurged,is.factor)==TRUE)
colnames(dataPurged)[-factor.variables] = paste0("accel.", colnames(dataPurged)[-factor.variables])

freq.1 = dataPurged
freq.2 = dataPurged
freq.3 = dataPurged

diff.1 = dataPurged
diff.2 = dataPurged
accel = dataPurged

freq.1 %>% gather(key = variable, value = value, -success,-name, -antigentype) -> freq.1.l
freq.3 %>% gather(key = variable, value = value, -success,-name, -antigentype) -> freq.3.l
freq.2 %>% gather(key = variable, value = value, -success,-name, -antigentype) -> freq.2.l
diff.2 %>% gather(key = variable, value = value, -success, -name, -antigentype) -> gp.2.l
diff.1 %>% gather(key = variable, value = value, -success, -name, -antigentype)  -> gp.1.l
accel %>% gather(key = variable, value = value, -success, -name, -antigentype) -> accel.l

the.whole.data = bind_rows(freq.1.l, freq.2.l, freq.3.l, gp.2.l, gp.1.l, accel.l) %>%
  spread(key = variable, value = value) %>% 
  select(-antigentype)

trial.multiple <- lme4::glmer(success ~.-name + (1 | name), data = the.whole.data, family = binomial,
                            control = glmerControl(optimizer="bobyqa"),nAGQ=0)

vif.results = data.frame(vif.mer(trial.multiple))
View(vif.results)

#1/22 
vif.excluded = c("freq.2.individual.meanMut", "freq.2.individual.meanBeta", "freq.1.individual.meanMut", "freq.3.individual.meanMut",
                 "freq.1.individual.meanBeta", "gp.1.individual.meanBeta", "gp.1.individual.meanMut", "gp.2.individual.meanBeta",
                 "accel.individual.meanMut", "freq.2.individual.varMut", "freq.2.individual.varBeta", "freq.2.individual.varR", "freq.3.individual.varMut",
                 "gp.1.individual.varMut", "accel.individual.varMut", "freq.2.meanBeta", "freq.1.individual.varBeta", "freq.1.individual.varR",
                 "gp.1.individual.varBeta", "gp.1.individual.varR","gp.2.individual.varBeta", "gp.2.individual.varR", "freq.2.ratio.meanBeta",
                 "freq.2.meanSigma", "freq.3.meanBeta", "freq.2.meanR", "freq.2.individual.meanR", "gp.2.individual.varMut", "freq.2.ratio.meanSigma",
                 "freq.2.individual.meanSigma", "freq.2.ratio.meanR", "freq.3.ratio.meanBeta", "freq.3.individual.meanBeta", "freq.3.meanSigma",
                 "freq.1.individual.meanR", "gp.1.individual.meanR","gp.2.individual.meanR", "freq.3.ratio.meanSigma", "gp.1.ratio.meanSigma",
                 "accel.ratio.meanSigma", "freq.3.meanR", "gp.1.meanR", "accel.ratio.meanBeta", "freq.1.individual.meanSigma", "accel.meanR",
                 "accel.individual.meanBeta", "freq.3.individual.meanSigma","accel.meanSigma", "accel.ratio.meanR", "freq.3.ratio.meanR", "gp.2.ratio.meanR",
                 "freq.2.covBetaSigma", "freq.2.varSigma", "freq.3.individual.meanR", "gp.2.ratio.meanBeta", "accel.individual.meanSigma", "freq.2.individual.varSigma",
                 "accel.individual.meanR", "freq.1.ratio.meanSigma", "gp.2.meanR", "freq.2.varR", "freq.1.covBetaSigma", "freq.1.meanSigma", 
                 "gp.2.ratio.meanSigma", "freq.2.meanLoad", "freq.1.varSigma", "gp.1.varSigma", "gp.2.varSigma", "gp.1.ratio.meanR", "freq.3.covBetaSigma",
                 "freq.3.individual.varSigma", "accel.varSigma","freq.2.diversity", "freq.2.antigenicDiversity", "freq.1.meanBeta", "freq.2.ratio.mutation",
                 "freq.3.meanLoad", "freq.2.ratio.varBeta", "freq.2.entropy", "freq.2.ratio.varR", "freq.3.diversity", "freq.3.antigenicDiversity",
                 "freq.2.antigenicTypes", "freq.2.dominant.freq", "freq.2.varBeta","freq.2.prop.I", "freq.2.totalI", "freq.2.totalS", "freq.1.individual.varMut",
                 "freq.1.ratio.mutation","accel.individual.varBeta", "accel.individual.varR", "freq.3.individual.varBeta", "freq.3.individual.varR",
                 "freq.3.entropy", "accel.ratio.mutation", "freq.1.ratio.varBeta", "gp.1.ratio.varBeta", "gp.2.ratio.varBeta", "freq.3.antigenicTypes",
                 "freq.1.antigenicTypes", "freq.3.dominant.freq", "gp.1.antigenicDiversity", "freq.1.varBeta", "freq.3.prop.I", "freq.3.totalI", "freq.3.totalS",
                 "gp.1.prop.I", "gp.1.totalI", "gp.1.totalS", "gp.1.meanLoad", "gp.1.diversity", "freq.2.ratio.varSigma", "freq.3.ratio.varR","gp.1.ratio.varR",
                 "accel.ratio.varR", "freq.1.ratio.meanBeta", "gp.1.antigenicTypes", "freq.2.tmrca", "gp.2.ratio.mutation", "freq.3.varSigma", "gp.1.ratio.mutation",
                 "gp.1.entropy", "gp.1.dominant.freq", "freq.1.entropy", "gp.2.ratio.varR", "accel.antigenicDiversity", "gp.2.antigenicTypes", "gp.1.meanBeta")

the.whole.data.purged = the.whole.data[,-which(colnames(the.whole.data) %in% vif.excluded)]
trial.purged <- lme4::glmer(success ~.-name + (1 | name), data = the.whole.data.purged, family = binomial,
                            control = glmerControl(optimizer="bobyqa"),nAGQ=0)
vif.results = data.frame(vif.mer(trial.purged))
View(vif.results)

dataPurged = the.whole.data.purged
write.csv(dataPurged, file = "../results/dataPurged.0129.csv")





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





  
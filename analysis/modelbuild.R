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



pop.level = c("day", "entropy","infected", "diversity", "antigenicDiversity", "antigenicTypes",
              "totalI", "meanR", "tmrca", "totalS", "dominant.freq")
#########################################################################
freq.df.subset %>%
  filter(freq == "freq.02") -> freq.2.subset

data = freq.3.subset

data  <- within(data, {
  success <- factor(success) #, levels=c("Est.", "Transient"))
  name <- as.factor(name)
  freq <- as.factor(freq)
  antigentype <- as.factor(antigentype)
})
 
#dataScaled = data %>% select(-simDay, -netau, -infected) %>%
#  mutate_if(is.numeric, scale)# Snap shot 
dataScaled = data %>% select(success, name, antigentype, one_of(pop.level), -day) %>%
  mutate_if(is.numeric, scale)# Snap shot 
dataPurged = dataScaled

########## Fit first model 
trial.all <- glmer(success ~.-name-antigentype + (1| name), data = dataScaled, family = binomial,
                   control = glmerControl(optimizer="bobyqa"),nAGQ=0)
vif.results = data.frame(vif.mer(trial.all))
View(vif.results)

# Snap Shot 
vif.excluded = c("antigenicTypes", "dominant.freq")
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
colnames(dataPurged)[-factor.variables] = paste0("freq1.", colnames(dataPurged)[-factor.variables])

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
vif.excluded = c("freq2.diversity", "freq2.antigenicDiversity", "freq2.entropy", "freq2.antigenicTypes",
                 "freq3.diversity", "freq3.antigenicDiversity", "freq2.infected", "freq2.dominant.freq",
                 "freq3.entropy", "freq2.meanR", "freq3.antigenicTypes", "freq1.antigenicTypes", "freq2.totalI",
                 "freq2.totalS", "freq3.dominant.freq", "freq2.tmrca", "freq3.infected", "gp1.infected",
                 "gp1.antigenicTypes", "accel.antigenicDiversity", "gp1.diversity", "gp1.dominant.freq",
                 "freq1.dominant.freq", "gp1.entropy", "freq1.meanR", "gp1.meanR", "gp2.antigenicTypes")


the.whole.data.purged = the.whole.data[,-which(colnames(the.whole.data) %in% vif.excluded)]
trial.purged <- lme4::glmer(success ~.-name + (1 | name), data = the.whole.data.purged, family = binomial,
                            control = glmerControl(optimizer="bobyqa"),nAGQ=0)
vif.results = data.frame(vif.mer(trial.purged))
View(vif.results)

dataPurged = the.whole.data.purged
write.csv(dataPurged, file = "../results/dataPurged.0207.realword.csv")





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





  
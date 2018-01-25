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


#########################################################################################
freq.df.subset %>%
  filter(freq == "freq.03") -> freq.3.subset

data = accelerate

data  <- within(data, {
  success <- factor(success) #, levels=c("Est.", "Transient"))
  name <- as.factor(name)
#  freq <- as.factor(freq)
  antigentype <- as.factor(antigentype)
})

# changing here -- used to take out totalI 
dataScaled = data %>% select(-simDay, -netau, -infected) %>%
  mutate_if(is.numeric, scale)# Snap shot 
dataPurged = dataScaled

########## Fit first model 
trial.all <- glmer(success ~.-name-antigentype + (1| name), data = dataScaled, family = binomial,
                   control = glmerControl(optimizer="bobyqa"),nAGQ=0)
vif.results = data.frame(vif.mer(trial.all))
View(vif.results)

# Snap Shot 
#vif.excluded = c("individual.meanMut", "individual.varMut", "individual.meanBeta", "individual.meanR", "individual.meanSigma",
#                 "ratio.meanSigma", "meanSigma", "covBetaSigma", "meanBeta", "individual.varR", "individual.varBeta", "antigenicTypes",
#                 "ratio.varR")

# 1/22/
vif.excluded = c("individual.meanMut", "individual.varMut", "individual.meanBeta", "individual.meanR", "individual.meanSigma",
                 "ratio.meanSigma", "meanSigma", "covBetaSigma", "meanBeta", "individual.varR", "individual.varBeta", "ratio.varR")

# 1/22/ Growth 
vif.excluded = c("individual.meanMut", "individual.varBeta", "individual.varR", "ratio.meanR", "meanR", "ratio.meanBeta",
                 "individual.meanR", "ratio.meanSigma", "varSigma", "individual.varMut")

dataPurged = dataScaled[,-which(colnames(dataScaled) %in% vif.excluded)]
trial.purged <- lme4::glmer(success ~.-name-antigentype + (1 | name), data = dataPurged, family = binomial,
                   control = glmerControl(optimizer="bobyqa"),nAGQ=0)
vif.results = data.frame(vif.mer(trial.purged))
View(vif.results)
vif.results %>% arrange(desc(vif.mer.trial.purged.)) %>% head()
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
#vif.excluded = c("freq.2.ratio.meanBeta", "freq.1.ratio.meanBeta", "freq.2.meanR", "diff.1.meanSigma",
#                 "freq.1.meanR", "diff.1.individual.meanR", "diff.1.individual.meanSigma", "diff.2.individual.meanSigma",
#                 "freq.2.varSigma", "freq.2.ratio.meanR", "freq.3.varSigma", "diff.1.ratio.varR", "accel.ratio.mutation",
#                 "accel.ratio.varR", "diff.1.individual.meanBeta", "accel.ratio.varBeta", "accel.individual.meanBeta",
#                 "freq.1.ratio.mutation", "freq.2.ratio.mutation", "accel.meanLoad", "freq.3.ratio.mutation")

vif.excluded = c("freq.1.individual.meanMut", "freq.1.individual.meanBeta", "freq.2.individual.meanMut", "freq.3.individual.meanMut", 
                "diff.1.individual.meanMut", "diff.1.individual.meanBeta", "diff.2.individual.meanMut", "accel.individual.meanMut", 
                "diff.1.individual.varMut", "diff.1.individual.varBeta", "diff.1.individual.varR", "accel.individual.varMut",
                "accel.individual.varBeta", "accel.individual.varR", "freq.2.individual.varMut", "freq.3.individual.varMut",
                "freq.1.individual.varMut", "diff.2.individual.varMut", "diff.1.ratio.meanR", "accel.ratio.meanR", "freq.2.individual.meanBeta",
                "freq.3.individual.meanBeta", "diff.1.individual.meanR", "diff.1.ratio.meanSigma", "accel.ratio.meanSigma", "freq.2.ratio.meanSigma",
                "freq.3.ratio.meanSigma", "diff.1.meanR", "diff.2.individual.meanR", "freq.2.individual.meanR", "freq.1.individual.meanR",
                "accel.ratio.meanBeta", "accel.meanR", "freq.2.meanR", "freq.3.meanR", "freq.1.individual.meanSigma", "freq.2.individual.meanSigma",
                "freq.3.individual.meanSigma", "diff.2.ratio.meanR", "freq.2.ratio.meanR", "accel.individual.meanBeta", "accel.meanSigma",
                "freq.3.individual.meanR", "diff.2.ratio.meanBeta", "freq.1.ratio.meanBeta", "freq.3.ratio.meanBeta", "accel.individual.meanSigma",
                "freq.1.ratio.meanR", "diff.2.meanR", "diff.2.ratio.meanSigma", "freq.1.meanSigma", "freq.2.meanSigma", "diff.1.varSigma",
                "freq.3.meanSigma", "freq.1.ratio.meanSigma", "accel.varSigma", "freq.2.varSigma", "freq.3.varSigma", "accel.individual.meanR",
                "accel.varR", "freq.1.covBetaSigma", "freq.2.covBetaSigma", "freq.3.covBetaSigma", "diff.2.varSigma", "freq.1.meanBeta",
                "freq.2.meanBeta", "freq.3.meanBeta", "accel.ratio.varBeta", "freq.2.individual.varBeta", "freq.2.individual.varR", 
                "freq.3.individual.varBeta", "freq.3.individual.varR", "freq.1.individual.varBeta", "freq.1.individual.varR", "diff.1.ratio.varR",
                "diff.2.individual.varBeta", "diff.2.individual.varR", "accel.ratio.varR", "accel.ratio.mutation", "freq.2.ratio.varR",
                "freq.3.ratio.varR", "freq.1.ratio.varR", "diff.1.meanLoad", "freq.1.ratio.mutation", "freq.2.ratio.meanBeta", "diff.2.ratio.mutation",
                "freq.2.ratio.mutation", "diff.2.ratio.varR", "diff.1.ratio.mutation", "freq.1.entropy", "diff.1.ratio.varSigma",
                "freq.2.entropy", "freq.3.entropy", "freq.1.varSigma", "accel.antigenicDiversity", "accel.diversity", "diff.2.ratio.varSigma", 
                "freq.2.ratio.varSigma", "diff.1.entropy", "diff.1.dominant.freq", "diff.1.antigenicTypes")


the.whole.data.purged = the.whole.data[,-which(colnames(the.whole.data) %in% vif.excluded)]
trial.purged <- lme4::glmer(success ~.-name + (1 | name), data = the.whole.data.purged, family = binomial,
                            control = glmerControl(optimizer="bobyqa"),nAGQ=0)
vif.results = data.frame(vif.mer(trial.purged))
View(vif.results)

dataPurged = the.whole.data.purged
write.csv(dataPurged, file = "../results/lower.vif.prop.csv")





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





  
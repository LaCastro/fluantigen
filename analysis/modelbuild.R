library(plyr)
library(dplyr)
library(tidyverse)
library(broom)
library(lme4)
library(pROC)
library(caret)
library(cvTools)
library(cowplot)
library(data.table)
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

data = freq.df.subset

data  <- within(data, {
  success <- factor(success) #, levels=c("Est.", "Transient"))
  name <- as.factor(name)
  freq <- as.factor(freq)
  antigentype <- as.factor(antigentype)
})
 
dataScaled = data %>% select(-simDay, -netau, -infected, -day, -freq) %>%
  mutate_if(is.numeric, scale) # Snap shot 
#dataScaled = data %>% select(success, name,  -day) %>%
#  mutate_if(is.numeric, scale)# Snap shot 

########## Fit first model 

select_variables <- function(data.set, data.names) {
  colinear = TRUE
  vif.excluded = vector()
  while(colinear == TRUE) { 
    trial.all <- glmer(success ~.-name-antigentype + (1| name), data = data.set, family = binomial,
                       control = glmerControl(optimizer="bobyqa"),nAGQ=0)
    vif.results = data.frame(vif.mer(trial.all))
    setDT(vif.results, keep.rownames = TRUE)
    candidate.row = vif.results %>% 
      arrange(desc(vif.mer.trial.all.)) %>% 
      slice(1)
    
    if(candidate.row$vif.mer.trial.all. > 5) {
      vif.excluded = c(vif.excluded, candidate.row$rn)
      data.set =  data.set[,-which(colnames(data.set) %in% vif.excluded)]
    } else {
      colinear = FALSE
    }
  }
  return(list(data = data.set, excluded.variables = vif.excluded))
}
eliminated.vif = select_variables(dataScaled, colnames(dataScaled))

eliminated.variables = eliminated.vif[[2]]
dataPurged = eliminated.vif[[1]] %>% select(-antigentype)

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





  
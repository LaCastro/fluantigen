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
pop.level = c("day", "entropy","infected", "diversity", "antigenicDiversity", "antigenicTypes",
              "totalI", "meanR", "tmrca", "totalS", "dominant.freq")
#########################################################################
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
select_variables_not_ant <- function(data.set, data.names) {
  colinear = TRUE
  vif.excluded = vector()
  while(colinear == TRUE) { 
    trial.all <- glmer(success ~.-name + (1| name), data = data.set, family = binomial,
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
#########################################################################

### CHANGE IF SNAPSHOT 
desired.freq = .01
type = "rf"
##############

freq.df.subset %>%
  filter(freq == desired.freq) -> freq.at.rf 
data = freq.at.rf 
if(type == "rf") {
  data  <- within(data, {
    success <- factor(success) #, levels=c("Est.", "Transient"))
    name <- as.factor(name)
    freq <- as.factor(freq)
    antigentype <- as.factor(antigentype)
  })
  dataScaled = data %>% select(-simDay, -netau, -infected) %>%
    mutate_if(is.numeric, scale) # Snap shot 
  
}
if(type == "egf") {
  data  <- within(data, {
    success <- factor(success) #, levels=c("Est.", "Transient"))
    name <- as.factor(name)
    antigentype <- as.factor(antigentype)
  })
  dataScaled = data %>% select(-simDay, -netau, -infected, -freq, -day) %>%
    mutate_if(is.numeric, scale) # Snap shot 
}

########## Fit first model 
#eliminated.vif = select_variables(dataScaled, colnames(dataScaled))
#eliminated.variables = eliminated.vif[[2]]
#dataPurged = eliminated.vif[[1]] %>% select(-antigentype)

################# Combing multiple time point data sets 
factor.variables = which(sapply(dataScaled,is.factor)==TRUE)
colnames(dataScaled)[-factor.variables] = paste0("freq", desired.freq, colnames(dataScaled)[-factor.variables])

# Snapshot 
freq.1 = dataScaled
freq.2 = dataScaled
freq.3 = dataScaled

diff.1 = dataScaled
diff.2 = dataScaled

accelerate.1 = dataScaled


############ Making it tidy
freq.1 %>% gather(key = variable, value = value, -success,-name, -antigentype) -> freq.1.l
freq.2 %>% gather(key = variable, value = value, -success,-name, -antigentype) -> freq.2.l
freq.3 %>% gather(key = variable, value = value, -success,-name, -antigentype) -> freq.3.l

diff.2 %>% gather(key = variable, value = value, -success, -name, -antigentype) -> gp.2.l
diff.1 %>% gather(key = variable, value = value, -success, -name, -antigentype)  -> gp.1.l

accelerate.1 %>% gather(key = variable, value = value, -success, -name, -antigentype) -> accel.1.l

the.whole.data = bind_rows(freq.1.l, freq.2.l, freq.3.l, gp.2.l, gp.1.l, accel.1.l) %>%
  spread(key = variable, value = value) %>% 
  select(-antigentype)


eliminated.vif = select_variables_not_ant(the.whole.data, colnames(the.whole.data))
eliminated.variables = eliminated.vif[[2]]
dataPurged = eliminated.vif[[1]]

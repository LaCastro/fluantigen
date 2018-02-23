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
library(magrittr)
#################################################################
########################################################################
select_variables <- function(data.set) {
  colinear = TRUE
  vif.excluded = vector()
  while(colinear == TRUE) { 
    trial.all <- glmer(success ~.-trial-antigentype-group + (1| trial), data = data.set, family = binomial,
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
select_variables_not_ant <- function(data.set) {
  colinear = TRUE
  vif.excluded = vector()
  while(colinear == TRUE) { 
    trial.all <- glmer(success ~.-trial + (1| trial), data = data.set, family = binomial,
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
desired.time = "diff"

data = full.data.df %>% filter(time == desired.time)

data  <- within(data, {
    success <- factor(success) #, levels=c("Est.", "Transient"))
    trial <- as.factor(trial)
    group <- as.factor(group)
    antigentype <- as.factor(antigentype)
    variable <- as.factor(variable)
    value <- as.numeric(as.character(value))
  })

data %>% 
  select(antigentype, success, variable, value, trial, group) %>%
  spread(key = variable, value = value) %>% select(-netau) %>%
  mutate(id = paste0(antigentype, "_",group)) %>% 
  mutate_at("id", as.factor) %>%
  mutate_if(is.numeric, scale) -> dataScaled

########## Fit first model 
#eliminated.vif = select_variables(dataScaled)
#eliminated.variables = eliminated.vif[[2]]
#dataPurged = eliminated.vif[[1]] %>% select(-antigentype)

############# Combing multiple time point data sets 
factor.variables = which(sapply(dataScaled,is.factor)==TRUE)
colnames(dataScaled)[-factor.variables] = paste0(desired.time, colnames(dataScaled)[-factor.variables])

# Snapshot 
time.1 = dataScaled
time.2 = dataScaled 
diff.1 = dataScaled

############ Making it tidy
time.1 %>% select(-group, -antigentype,-t1day) %>% gather(key = variable, value = value, -success, -trial, -id) -> time.1.l
time.2 %>% select(-group, -antigentype,-t2day) %>% gather(key = variable, value = value, -success, -trial, -id) -> time.2.l
diff.1 %>% select(-group, -antigentype) %>% gather(key = variable, value = value, -success, -trial, -id) -> diff.1.l
the.whole.data = bind_rows(time.1.l, time.2.l, diff.1.l) %>%
  spread(key = variable, value = value) %>% 
  select(-id)

eliminated.vif = select_variables_not_ant(the.whole.data)
eliminated.variables = eliminated.vif[[2]]
dataPurged = eliminated.vif[[1]]

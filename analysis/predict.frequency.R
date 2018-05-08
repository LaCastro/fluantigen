rm(list=ls())
source('loadLibraries.R')
source('analysis_functions.R')
source('plotting_functions.R')
library(afex)
library(Metrics)
library(MuMIn)


select_variables_frequency <- function(data.set, time.ahead) {

  colinear = TRUE
  vif.excluded = vector()
  excluded.time = time.points[-which(time.points == time.ahead)]
  data.sub = data.set %>% select(-antigentype, -success, -group) %>% 
    select(-one_of(excluded.time))
  
  while(colinear == TRUE) { 
    
    formula = paste0(time.ahead, "~.-trial + (1 | trial)")
    trial.all <- lmer(formula, data = data.sub) # Change 
    vif.results = data.frame(vif.mer(trial.all))
    setDT(vif.results, keep.rownames = TRUE)
    candidate.row = vif.results %>% 
      arrange(desc(vif.mer.trial.all.)) %>% 
      slice(1)
    
    if(candidate.row$vif.mer.trial.all. > 5) {
      vif.excluded = c(vif.excluded, candidate.row$rn)
      data.sub =  data.sub[,-which(colnames(data.sub) %in% vif.excluded)]
    } else {
      colinear = FALSE
    }
  }
  return(list(data = data.sub, excluded.variables = c(vif.excluded,excluded.time)))
}
get_final_model <- function(results, data , dependent.variable, part) {
  formula = paste0(dependent.variable,"~", paste0(results, collapse = "+"), fixed.part.2)
  if(part == "part1") { 
    model.null = lmer_alt(formula, family = binomial(link = "logit"), data = data)
  }
  else {
    data %<>% filter(above.zero == 1)
    model.null = lmer_alt(formula, data = data)
  }
  return(model.null)
}
## Func
calculate.day.ahead = function(sample.time, time.ahead) {
  sample.t = sample.time + time.ahead
  return(sample.t)
}
get_frequencies_time_ahead <- function(antigen.frequency.df,sample.times.l) {
  
  days = sample.times.l$day
  antigen.frequency.df %>%
    filter(day %in% days) %>%
    rename(day.ahead = day) %>%
    select(-infected) %>%
    left_join(sample.times.l, by = c("day.ahead" = "day"))  %>%
    mutate_at("antigentype", as.character)
} 
join_timeahead = function(surveillance.data, frequency) { 
  right_join(surveillance.data, frequency, by = c("antigentype", "group"))
}

############################################################################


# Step 1. Put together list trials of transient and successful antigens ; this will be based on max.freq, and days above
tropics.folder = "../data/tropics/eligible/"
trial.dirs = dir(tropics.folder)

# This is a slow step
antigen.frequencies = map(trial.dirs, function(x) read.table(paste0(tropics.folder,x,"/out.antigenFrequencies.txt"), header = TRUE)) 
names(antigen.frequencies) = trial.dirs

######## Step 2: For each list entry need to calculate max frequency and days above
thres = .2

# Calculating the two criteria for determining success and if successful how long was it there 
antigen.frequencies %>%
  map(function(x) calculate_days_above(x, threshold = .2)) -> days.above

antigen.frequencies %>%
  map(find_max_frequency) %>%
  map(function(x) mutate_at(x, "antigentype", as.character)) -> max.frequencies

# Creating the combined data frame and using the criteria to assign a label 
rel.frequency.thres = .01 # Threshold at which we start monitoring/expect to pick up in surveillance 
full.data = map2(max.frequencies, days.above, left_join) %>%
  map(replace_na_zeros) %>% 
  map(determine_success_labels, max.rf = rel.frequency.thres)

####### Step 3: Subsetting for just transient and yes; and getting rid of zeros 
subset.data = map(full.data, filter_out_loss) 
subset.df = do.call("rbind", subset.data)
subset.df$name = rep(trial.dirs, sapply(subset.data, nrow))

subset.df %>%
  filter(antigentype != 0) -> subset.analyze

######### Step 4 - Calculate Entropy
entropy = map(antigen.frequencies, calculate_entropy)
names(entropy) = trial.dirs

######### Step 5 Generate Sample Times
set.seed(103114)
possible.sample = seq(from = 1833, to = 10947, by = 14)
sample.time.one = sample(possible.sample, 10, replace=FALSE) 
sample.time.two = adply(sample.time.one, .margins = 1, generate_second_time,.id = NULL)
sample.times = data.frame(t1 = sample.time.one, t2= sample.time.two$V1) %>% arrange(t1)

######################## Function to collect population level at those times
# For trial time 1 
data.1.list = map(trial.dirs, create_trial_dataset, sample.times$t1) 
data.1.list = map(data.1.list, create_time_variable)

sample.times %>%
  mutate(three.months = calculate.day.ahead(t1, 6*14),
         six.months = calculate.day.ahead(t1, 12*14),
         nine.months = calculate.day.ahead(t1, 18*14),
         twelve.months = calculate.day.ahead(t1, 24*14)) -> sample.times

sample.times %>%
  mutate(group = row_number()) %>%
  gather(key = time, value = day, -t1, -t2, -group) %>%
  select(-t1, -t2) -> sample.times.l

frequency.time.ahead = map(antigen.frequencies, get_frequencies_time_ahead, sample.times.l)
trial.combine = map2(data.1.list, frequency.time.ahead, join_timeahead) %>% 
  map(function(x) x %>% mutate(time.stamp = 't1')) %>%
  dplyr::bind_rows(.id = 'list')

trial.combine %>% 
  filter(!is.na(success)) %>%
  select(-day.ahead) %>%
  spread(key = time, value = frequency, fill = 0) -> trial.combine.w

remove.variables = c("date", "day", "day.success", "stage", "simDay")
trial.combine.w %<>% filter(!(variable %in% remove.variables))

########################################################## Eliminating Variables 
data = trial.combine.w
data  <- within(data, {
  success <- factor(success) #, levels=c("Est.", "Transient"))
  trial <- as.factor(trial)
  group <- as.factor(group)
  antigentype <- as.factor(antigentype)
  variable <- as.factor(variable)
  value <- as.numeric(as.character(value))
})

time.points = c("three.months", "six.months", "nine.months", "twelve.months")

 data %>%
   select(antigentype, success, variable, value, trial, group, nine.months, six.months, twelve.months, three.months) %>%
   spread(key = variable, value = value) %>%
   select(-netau) %>%
   ggplot(aes(frequency)) + geom_histogram(binwidth = .001) + facet_wrap(~group)


data %>% 
  select(antigentype, success, variable, value, trial, group, nine.months, six.months, twelve.months, three.months) %>%
  spread(key = variable, value = value) %>% 
  select(-netau) %>%
  mutate(three.fold = three.months/frequency,
         six.fold = six.months/frequency,
         nine.fold = nine.months/frequency,
         twelve.fold = twelve.months/frequency) %>%
  mutate_at(-c(1:8, 41:44), scale) ->  dataScaled
  
# dataScaled %>%
#   select(nine.months, six.months, twelve.months, three.months) %>%
#   gather(key = time, value = value) %>%
#   ggplot(aes(value)) + geom_histogram() +facet_wrap(~time) -> histogram.frequency.advance
# save_plot(plot = histogram.frequency.advance,"exploratory.figs/histogram.frequency.advance.pdf", base_height = 8)

############ Chooose Which Month Interested In 
time.points = c("three.months", "six.months", "nine.months", "twelve.months", "three.fold", "six.fold", "nine.fold", "twelve.fold")
eliminated.vif = select_variables_frequency(dataScaled, time.ahead = "twelve.fold")
eliminated.variables = eliminated.vif[[2]]
dataPurged = eliminated.vif[[1]] 

######################################################### Modeling Building 
set.seed(62417)
folds.length = 5;
dataPurged %<>% rename(name = trial) 

folds = cvFolds(n = length(unique(dataPurged$name)), K = folds.length, type = "random")
#######################################################################################
#### Part One
# Fit model for above 0 or below
dependent.variable = "twelve.fold" ## Change here
dataPurged %>%
  mutate(above.zero = ifelse(twelve.fold > 0.0, 1,0)) %>% ### Change here
  mutate_at("above.zero", as.factor) -> dataPurged.2

build = TRUE
variable.list = c()
variable.aic = c()
variable.perf = c()

fixed.part.1 =  "above.zero ~ value  + "
fixed.part.2 = " + (1 | name)"

aic.values = rep(0, folds.length)
anova.results = rep(0,folds.length)
glmerperf = rep(0, folds.length)

while(build == TRUE) {
  # Step 1: Sets the formula and manipulates the data set
  if(length(variable.list) == 0 )  { # This is the first round
    formula = "above.zero ~ value + (1|name)"
    dataPurged.2 %>%
      ungroup() %>%
      gather(key = variable, value = value, -one_of(dependent.variable), -above.zero, -name) %>%
      mutate_at("value", as.numeric) -> data.scaled.l
  } else
    { # There are already significant variables
    variable.part = paste(variable.list, collapse = "+")
    formula = paste0(fixed.part.1, variable.part, fixed.part.2)
    formula.null = paste0("above.zero ~", paste0(variable.list, collapse = "+"), fixed.part.2)

    dataPurged.2 %>%
      ungroup() %>%
      gather(key = variable, value = value, -one_of(dependent.variable), -name, -above.zero, -one_of(variable.list)) %>%
      mutate_at("value", as.numeric) -> data.scaled.l
  }

  model.results = ddply(data.scaled.l, .var = c("variable"), .fun = function(x) {
    # K-fold validation for each model based on a single term
    for(n in 1:folds.length) {
      test.ids = folds$subsets[folds$which==n]
      test.trials = unique(data.scaled.l$name)[test.ids]
      testdata = x[which(x$name %in% test.trials), ]
      traindata = x[-which(x$name %in% test.trials),]
      model.new <- glmer(formula, data = traindata, family = binomial(link = 'logit'))

      aic.values[n] = glance(model.new)$AIC

      glmer.probs <- predict(model.new, newdata=testdata, type="response", allow.new.levels=TRUE)
      if (testdata %>% filter(above.zero == 0) %>% nrow() == 0) {
        glmerperf[n] = NA
      } else {
      glmer.ROC <- roc(predictor=glmer.probs, response=testdata$above.zero)
      glmerperf[n] <- glmer.ROC$auc
      }

      if (length(variable.list) > 0) { # Compare new model against the null
        model.null = lme4::glmer(formula.null, family = binomial(link = "logit"), data = traindata)
        anova.test = anova(model.null, model.new, refit=FALSE)
        anova.results[n] = anova.test$`Pr(>Chisq)`[2]
      }
    }

    aic.mean = mean(aic.values)
    perf.mean = mean(glmerperf, na.rm = TRUE)
    anova.mean = mean(anova.results)

    return(t(c(aic = aic.mean, perf = perf.mean, anova = anova.mean)))# put the se of the performance
  })

  variable = model.results %>% arrange(aic) %>% slice(1)  # Select the variable that is the best

  if(length(variable.list) == 0 ) { # Is this the first? If yes, add this variable
    variable.list = c(variable.list, variable$variable)
    variable.aic = c(variable.aic, variable$aic)
    variable.perf = c(variable.perf, variable$perf)
    print(paste0("Adding : ", variable$variable))
  } else if (variable$anova < .05) { # Determine if this new variable beats the null
    variable.list = c(variable.list, variable$variable)
    variable.aic = c(variable.aic, variable$aic)
    variable.perf = c(variable.perf, variable$perf)
    print(paste0("Adding : ", variable$variable))
  } else {
    print(paste0("Final Part 1 Model is: ", formula.null))
    build = FALSE
  }
}

first.part.results = bind_cols(variable = variable.list, per =  variable.perf, aic = variable.aic)
write.csv(first.part.results, '../results/050218.twelvemonths.part1.csv', row.names = FALSE)

part1.12.optimal.model = get_final_model(first.part.results$variable, data = dataPurged.2, dependent.variable = "above.zero", part = "part1")
saveRDS(part1.12.optimal.model, "../results/models/part1.threefold.optimal.months.rds")

###################### Second part
dataPurged.2 %>% filter(above.zero == 1) %>%
  mutate(twelve.fold.t = (twelve.fold)^(1/3)) %>%
  select(-twelve.fold) -> dataPurged.t
ggplot(dataPurged.t,aes(twelve.fold.t)) + geom_histogram()
dependent.variable = "twelve.fold.t"

build = TRUE
part2.variable.list = c()
part2.variable.aic = c()
part2.variable.rmse= c()

part2.fixed.part.1 =  paste0(dependent.variable, " ~ value  + ")
part2.fixed.part.2 = " + (1 | name)"

aic.values = rep(0, folds.length)
anova.results = rep(0,folds.length)
rmse.results = rep(0, folds.length)

while(build == TRUE) {
  # Step 1: Sets the formula and manipulates the data set
  if(length(part2.variable.list) == 0 )  { # This is the first round
    formula = paste(dependent.variable, "~value + (1 | name)")
    dataPurged.t %>%
      ungroup() %>%
      gather(key = variable, value = value, -one_of(dependent.variable), -above.zero, -name) %>%
      mutate_at("value", as.numeric) -> data.scaled.l
  } else
    { # There are already significant variables
    variable.part = paste(part2.variable.list, collapse = "+")
    formula = paste0(part2.fixed.part.1, variable.part, part2.fixed.part.2)
    formula.null = paste0(dependent.variable, "~", paste0(part2.variable.list, collapse = "+"), part2.fixed.part.2)

    dataPurged.t %>%
      ungroup() %>%
      gather(key = variable, value = value, -one_of(dependent.variable), -name, -above.zero, -one_of(part2.variable.list)) %>%
      mutate_at("value", as.numeric) -> data.scaled.l
  }

  model.results = ddply(data.scaled.l, .var = c("variable"), .fun = function(x) {
    # K-fold validation for each model based on a single term
    for(n in 1:folds.length) {
      test.ids = folds$subsets[folds$which==n]
      test.trials = unique(data.scaled.l$name)[test.ids]
      testdata = x[which(x$name %in% test.trials), ]
      traindata = x[-which(x$name %in% test.trials),]

      model.linear <- lme4::lmer(formula, data = traindata, subset = above.zero == 1)
      aic.values[n] = glance(model.linear)$AIC

      testdata %>% filter(above.zero == 1) -> test.sub
      lmer.probs <- predict(model.linear, newdata=test.sub,  type="response", allow.new.levels=TRUE)
      rmse.results[n] <- rmse(actual = test.sub[,dependent.variable], predicted = lmer.probs)
      if (length(part2.variable.list) > 0) { # Compare new model against the null
        model.null = lme4::lmer(formula.null, data = traindata, subset = above.zero == 1)
        anova.test = anova(model.null, model.linear, refit=FALSE)
        anova.results[n] = anova.test$`Pr(>Chisq)`[2]
      }
    }
    anova.mean = mean(anova.results)
    rmse.mean=mean(rmse.results)
    aic.mean = mean(aic.values)

    return(t(c(anova = anova.mean, rmse.mean = rmse.mean, aic = aic.mean)))# put the se of the performance
  })

  variable = model.results %>% arrange(aic) %>% slice(1)  # Select the variable that is the best

  if(length(part2.variable.list) == 0 ) { # Is this the first? If yes, add this variable
    part2.variable.list = c(part2.variable.list, variable$variable)
    part2.variable.aic = c(part2.variable.aic, variable$aic)
    part2.variable.rmse = c(part2.variable.rmse, variable$rmse.mean)

    print(paste0("Adding : ", variable$variable))
  } else if (variable$anova < .05) { # Determine if this new variable beats the null
    part2.variable.list = c(part2.variable.list, variable$variable)
    part2.variable.aic = c(part2.variable.aic, variable$aic)
    part2.variable.rmse = c(part2.variable.rmse, variable$rmse.mean)
    print(paste0("Adding : ", variable$variable))
  } else {
    print(paste0("Final Part 2 Model is: ", formula.null))
    build = FALSE
  }
}
part2.results = bind_cols(variable = part2.variable.list, variable.aic =  part2.variable.aic,  variable.rmse =  part2.variable.rmse)
write.csv(part2.results, '../results/041018.threemonths.part2.optimal.cubedroot.csv', row.names = FALSE)

#variables = c("individual.meanR", "frequency", "meanR", "ratio.varBeta", "varR", "ratio.meanBeta")
part2.model = get_final_model(part2.results$variable, data = dataPurged.t, dependent.variable = dependent.variable, part = "part2") #part2.results$variable
saveRDS(part2.model, "../results/models/part2.twelve.fold.optimal.cubedroot.rds")

#### Diagnostics
plot(part2.model)
car::qqPlot(residuals(part2.model), main = "Twelve Months - Cubed Root")

# 
################################
# 
# ############################ Testing 
# 
# set.seed(22591)
# possible.sample = seq(from = 1833, to = 10947, by = 14)
# test.time.one = sample(possible.sample, 5, replace=FALSE) 
# test.time.two = adply(test.time.one, .margins = 1, generate_second_time,.id = NULL)
# test.times = data.frame(t1 = test.time.one, t2= test.time.two$V1) %>% arrange(t1)
# 
# ##### Collect population level at those times
# # For trial time 1 
# test.list = map(trial.dirs, create_trial_dataset, test.times$t1) 
# test.list = map(test.list, create_time_variable)
# # Only need to do this the first time 
# sample.times %>%
#   mutate(two.months = calculate.day.ahead(t1, 4*14),
#          six.months = calculate.day.ahead(t1, 12*14),
#          nine.months = calculate.day.ahead(t1, 18*14),
#          twelve.months = calculate.day.ahead(t1, 26*14)) -> sample.times
# 
# sample.times %>%
#   mutate(group = row_number()) %>%
#   gather(key = time, value = day, -t1, -t2, -group) %>%
#   select(-t1, -t2) -> sample.times.l
# 
# frequency.time.ahead = map(antigen.frequencies, get_frequencies_time_ahead, sample.times.l)
# test.combine = map2(test.list, frequency.time.ahead, join_timeahead) %>% 
#   map(function(x) x %>% mutate(time.stamp = 't1')) %>%
#   dplyr::bind_rows(.id = 'list')
# 
# test.combine %>% 
#   filter(!is.na(success)) %>%
#   select(-day.ahead) %>%
#   spread(key = time, value = frequency, fill = 0) -> test.combine.w
# 
# remove.variables = c("date", "day", "day.success", "stage", "simDay")
# test.combine.w %<>% filter(!(variable %in% remove.variables))
# 
# ########################################################## Eliminating Variables 
# validate.data = test.combine.w
# validate.data  <- within(validate.data, {
#   success <- factor(success) #, levels=c("Est.", "Transient"))
#   trial <- as.factor(trial)
#   group <- as.factor(group)
#   antigentype <- as.factor(antigentype)
#   variable <- as.factor(variable)
#   value <- as.numeric(as.character(value))
# })
# 
# time.points = c("two.months", "six.months", "nine.months", "twelve.months")
# 
# validate.data %>% 
#   select(antigentype, success, variable, value, trial, group, nine.months, six.months, twelve.months, two.months) %>%
#   spread(key = variable, value = value) %>% 
#   select(-netau) %>%
#   mutate_at(-c(1:8), scale) ->  validate.dataScaled
# 
# 
# ## With data, have to transform to the desired model 
# validate.dataScaled %>% 
#   rename(name = trial) %>%
#   mutate(above.zero = ifelse(twelve.months > 0, 1,0)) -> dataScaled.zero
# 
# # part 1 test 
# test.probs <- predict(part1.model, newdata=dataScaled.zero, type="response", allow.new.levels=TRUE)
# test.ROC <- roc(predictor=test.probs, response=dataScaled.zero$above.zero)
# test.ROC$auc
# roc.values = data.frame(cbind(sen = test.ROC$sensitivities,
#                               spec = test.ROC$specificities,
#                               thres = test.ROC$thresholds))
# write.csv(roc.values, file = "../results/roc.twelvemonths.cubedroot.csv", row.names = FALSE)
# output.part1 = cbind(pred = test.probs, actual = dataScaled.zero$above.zero, success = dataScaled.zero$success, antigen = dataScaled.zero$antigentype, name = dataScaled.zero$name)
# write.csv(output.part1, file = "../results/output.part1.twelvemonths.cubedroot.csv", row.names = FALSE)
# 
# # part 2 test 
# dataScaled.zero %>% filter(above.zero == 1) %>% 
#   mutate(twelve.months.t = (twelve.months)^(1/3)) -> dataScaled.t
# 
# lmer.probs <- predict(part2.model, newdata=dataScaled.t,  type="response", allow.new.levels=TRUE)
# lmerRmse <- rmse(actual = dataScaled.t$twelve.months.t, predicted = lmer.probs) 
# output.part2 = cbind(pred = lmer.probs, actual = dataScaled.t$twelve.months.t, success= dataScaled.t$success, antigen = dataScaled.t$antigentype, name = dataScaled.t$name)
# write.csv(output.part2, file = "../results/output.part2.twelvemonths.cubedroot.csv", row.names = FALSE)
###################### 


### Fitting models with just the top three 

#eliminated.vif = select_variables_frequency(dataScaled, time.ahead = "three.months")

top.variables = c("trial", "three.months","six.months", "nine.months", "twelve.months",
"frequency", "success", "ratio.meanR", "varR")
dataPurged = dataScaled[,top.variables] 



######################################################### Modeling Building 
set.seed(62417)
folds.length = 5;
dataPurged %<>% rename(name = trial) 

folds = cvFolds(n = length(unique(dataPurged$name)), K = folds.length, type = "random")

#######################################################################################

## Trying to do a second model with two fits 

#### Part One
# Fit model for above 0 or below 
dependent.variable = "twelve.months" ## Change here
dataPurged %>% 
  mutate(above.zero = ifelse(three.months > 0.0, 1,0)) %>% ### Change here 
  mutate_at("above.zero", as.factor) -> dataPurged.2

build = TRUE
variable.list = c()
variable.aic = c()
variable.perf = c()

formula = "above.zero ~ ratio.meanR + frequency + varR + (1 |name)" 
glmerperf = rep(0, folds.length)


# K-fold validation for each model based on a single term 
for(n in 1:folds.length) { 
  test.ids = folds$subsets[folds$which==n]
  test.trials = unique(dataPurged.2$name)[test.ids]
  testdata = dataPurged.2[which(dataPurged.2$name %in% test.trials), ]
  traindata = dataPurged.2[-which(dataPurged.2$name %in% test.trials),]
  model.new <- glmer(formula, data = traindata, family = binomial(link = 'logit'))
  glmer.probs <- predict(model.new, newdata=testdata, type="response", allow.new.levels=TRUE)
  
  if (testdata %>% filter(above.zero == 0) %>% nrow() == 0) {
    glmerperf[n] = NA
  } else { 
    glmer.ROC <- roc(predictor=glmer.probs, response=testdata$above.zero)
    glmerperf[n] <- glmer.ROC$auc
  }
}
    
perf.three = glmerperf
perf.six = glmerperf    
perf.nine = glmerperf
perf.twelve = glmerperf
mean(perf.twelve)

first.part.results = bind_cols(variable = variable.list, per =  variable.perf, aic = variable.aic)
write.csv(first.part.results, '../results/030218.ninemonths.part1.csv', row.names = FALSE)

part1.3.model = get_final_model(results  = c("ratio.meanR", "frequency", "varR"), data = dataPurged.2, dependent.variable = "above.zero", part = "part1")

saveRDS(part1.model, "../results/models/part1.twelvemonths.rds")

###################### Second part 
dataPurged %>%
  mutate(three.months.t = (three.months)^(1/3),
         six.months.t = (six.months)^(1/3),
         nine.months.t = (nine.months)^(1/3),
         twelve.months.t = (twelve.months) ^(1/3)) %>%
  mutate(above.zero = ifelse(twelve.months > 0.0, 1,0)) %>%
  filter(above.zero == 1) -> dataPurged.t

ggplot(dataPurged.t,aes(twelve.months.t)) + geom_histogram()
dependent.variable = "twelve.months.t"

formula = paste0(dependent.variable, "~ratio.meanR + frequency + varR + (1|name)") 
rmse.results = rep(0, folds.length)
for(n in 1:folds.length) { 
  test.ids = folds$subsets[folds$which==n]
  test.trials = unique(dataPurged.t$name)[test.ids]
  testdata = dataPurged.t[which(dataPurged.t$name %in% test.trials), ]
  traindata = dataPurged.t[-which(dataPurged.t$name %in% test.trials),]
  
  model.linear <- lme4::lmer(formula, data = traindata, subset = above.zero == 1)
  
  testdata %>% filter(above.zero == 1) -> test.sub
  lmer.probs <- predict(model.linear, newdata=test.sub,  type="response", allow.new.levels=TRUE)
  rmse.results[n] <- rmse(actual = test.sub[,dependent.variable], predicted = lmer.probs) 
}

three.months.rmse = rmse.results; mean(three.months.rmse)
six.months.rmse = rmse.results
nine.months.rmse = rmse.results
twelve.months.rmse = rmse.results
  

part2.results = bind_cols(variable = part2.variable.list, variable.aic =  part2.variable.aic,  variable.rmse =  part2.variable.rmse)
write.csv(part2.results, '../results/030118.twomonths.part2.cubedroot.csv', row.names = FALSE)

variables = c("ratio.meanR", "frequency", "varR")
part2.3.otpimal.model = get_final_model(results = variables, data = dataPurged.t, dependent.variable = dependent.variable, part = "part2") #part2.results$variable
saveRDS(part2.model, "../results/models/part2.twelvemonths.cubedroot.rds")

#### Diagnostics
plot(part2.model)
car::qqPlot(residuals(part2.model), main = "Twelve Months - Cubed Root")

############################ Testing 


calculate_practicalResults = function(thres, prediction.data) {
  prediction.data %>%
    mutate(predict = ifelse(test.probs < thres, "0", "1"),
           thres = thres)
}

calculate_outcomeTypes = function(predictedOutcome.x, type) {
  predictedOutcome.x %<>%
    mutate(predict = as.factor(predict),
           truth = as.factor(truth))
  cm = confusionMatrix(reference = predictedOutcome.x$truth, data = predictedOutcome.x$predict, positive = "1")
  if(type == "class") {
    return(cm$byClass)
  } else {
    return(cm$overall)
  }
}


set.seed(22591)
possible.sample = seq(from = 1833, to = 10947, by = 14)
test.time.one = sample(possible.sample, 5, replace=FALSE) 
test.time.two = adply(test.time.one, .margins = 1, generate_second_time,.id = NULL)
test.times = data.frame(t1 = test.time.one, t2= test.time.two$V1) %>% arrange(t1)

##### Collect population level at those times
# For trial time 1 
test.list = map(trial.dirs, create_trial_dataset, test.times$t1) 
test.list = map(test.list, create_time_variable)
# Only need to do this the first time 
test.times %>%
  mutate(three.months = calculate.day.ahead(t1, 6*14),
         six.months = calculate.day.ahead(t1, 12*14),
         nine.months = calculate.day.ahead(t1, 18*14),
         twelve.months = calculate.day.ahead(t1, 24*14))  -> test.times

test.times %>%
  mutate(group = row_number()) %>%
  gather(key = time, value = day, -t1, -t2, -group) %>%
  select(-t1, -t2) -> test.times.l

frequency.time.ahead = map(antigen.frequencies, get_frequencies_time_ahead, test.times.l)
test.combine = map2(test.list, frequency.time.ahead, join_timeahead) %>% 
  map(function(x) x %>% mutate(time.stamp = 't1')) %>%
  dplyr::bind_rows(.id = 'list')

test.combine %>% 
  filter(!is.na(success)) %>%
  select(-day.ahead) %>%
  spread(key = time, value = frequency, fill = 0) -> test.combine.w

remove.variables = c("date", "day", "day.success", "stage", "simDay")
test.combine.w %<>% filter(!(variable %in% remove.variables))

########################################################## Eliminating Variables 
validate.data = test.combine.w
validate.data  <- within(validate.data, {
  success <- factor(success) #, levels=c("Est.", "Transient"))
  trial <- as.factor(trial)
  group <- as.factor(group)
  antigentype <- as.factor(antigentype)
  variable <- as.factor(variable)
  value <- as.numeric(as.character(value))
})

time.points = c("three.months", "six.months", "nine.months", "twelve.months")

validate.data %>% 
  select(antigentype, success, variable, value, trial, group, nine.months, six.months, twelve.months, three.months) %>%
  spread(key = variable, value = value) %>% 
  select(-netau) %>%
  mutate_at(-c(1:8), scale) ->  validate.dataScaled


## With data, have to transform to the desired model 
validate.dataScaled %>% 
  rename(name = trial) %>%
  mutate(above.zero = ifelse(three.months > 0, 1,0)) -> dataScaled.zero

# part 1 test 
test.probs <- predict(part1.3.optimal.model, newdata=dataScaled.zero, type="response", allow.new.levels=TRUE)
test.ROC <- roc(predictor=test.probs, response=dataScaled.zero$above.zero)
test.ROC$auc

potential.thres = seq(.01, .99, .01)
compare.prediction = data.frame(test.probs, truth = as.factor(dataScaled.zero$above.zero), success = dataScaled.zero$success, current = dataScaled.zero$frequency)
predictedOutcomes = map(potential.thres, calculate_practicalResults, compare.prediction)
performance.statistics = map(predictedOutcomes, calculate_outcomeTypes, "class")

performance.part1.3 = data.frame(do.call("rbind", performance.statistics))
performance.part1.12$thres = potential.thres

performance.part1.3 %<>% mutate(month = "three")
performance.part1.6 %<>% mutate(month = "six")
performance.part1.9 %<>% mutate(month = "nine")
performance.part1.12 %<>% mutate(month = "twelve")


performance.part1.3 %>%
  bind_rows(performance.part1.6, performance.part1.9,performance.part1.12) -> performance.part1
performance.part1$month = factor(performance.part1$month, levels = c("three", "six", "nine", "twelve"),
                                 labels = c("3", "6", "9", "12"))
performance.part1 %>%
  ggplot(aes(x = Sensitivity, y = Pos.Pred.Value, color = thres, shape = month)) + geom_point(size  = 3) +
  scale_color_viridis(breaks = seq(0,1,.2), option = "A", end = .9) +
  labs(color = "Probability \n Threshold", shape = "Months Ahead", y = "Positive Predictive Value") -> part1.compare.plot
save_plot(part1.compare.plot, filename = "exploratory.figs/part1.compare.plot.pdf", base_aspect_ratio = 1.3)


write.csv(roc.values, file = "../results/roc.twelvemonths.cubedroot.csv", row.names = FALSE)
output.part1 = cbind(pred = test.probs, actual = dataScaled.zero$above.zero, success = dataScaled.zero$success, antigen = dataScaled.zero$antigentype, name = dataScaled.zero$name)
write.csv(output.part1, file = "../results/output.part1.twelvemonths.cubedroot.csv", row.names = FALSE)

# Part 2 test 
validate.dataScaled %>% 
  rename(name = trial) %>%
  mutate(above.zero = ifelse(twelve.months > 0, 1,0)) -> dataScaled.zero


dataScaled.zero %>% filter(above.zero == 1) %>% 
  mutate(twelve.months.t = (twelve.months)^(1/3)) -> dataScaled.t


part2.12.optimal = readRDS("~/Dropbox/current_fluantigen/results/models/part2.twelvemonths.cubedroot.rds")
lmer.probs <- predict(part2.12.optimal, newdata=dataScaled.t,  type="response", allow.new.levels=TRUE)
lmerRmse <- rmse(actual = dataScaled.t$twelve.months.t, predicted = lmer.probs)
output.part2.12.optimal = data.frame(month = "twelve", pred = as.numeric(lmer.probs),
                                    actual = as.numeric(dataScaled.t$twelve.months.t), 
                                   success= dataScaled.t$success,
                            antigen = dataScaled.t$antigentype, name = dataScaled.t$name)

write.csv(output.part2.3.optimal, file = "../results/output.part2.3.optimal.cubedroot.csv", row.names = FALSE)


output.part2.3.optimal %>%
  mutate(residual = actual-pred) %>%
  ggplot(aes(x = pred^3, y = actual^3, color = success)) + geom_point() + geom_abline(slope = 1, intercept = 0) +
  scale_color_manual(values = c("black", "grey")) + labs(x = expression("True Frequency"), 
                                                         y = expression("Predicted Frequency"),
                                                         color = "Antigen Fate")  

############ Model Performance Combined
validate.dataScaled %>% 
  rename(name = trial) %>%
  mutate(above.zero = ifelse(three.months > 0, 1,0)) -> dataScaled.zero

# part 1 test 
test.probs <- predict(part1.3.model, newdata=dataScaled.zero, type="response", allow.new.levels=TRUE)
test.ROC <- roc(predictor=test.probs, response=dataScaled.zero$above.zero)
test.ROC$auc

potential.thres = seq(.01, .99, .01)
compare.prediction = data.frame(test.probs, truth = as.factor(dataScaled.zero$above.zero), success = dataScaled.zero$success)
predictedOutcomes = map(potential.thres, calculate_practicalResults, compare.prediction)
performance.statistics = map(predictedOutcomes, calculate_outcomeTypes, "class")
performance.part1.12 = data.frame(do.call("rbind", performance.statistics))
performance.part1.12$thres = potential.thres

compare.prediction %>%
  mutate(above.thres = ifelse(test.probs > optimal.six.thres$thres, 1, 0)) %>%
  group_by(success, truth) %>%
  count(above.thres) 
  
potential.figure5= plot_grid(part1.compare.plot, twelve.months.part2.plot)
save_plot(potential.figure5, filename = "exploratory.figs/potential.figure5.pdf", base_aspect_ratio = 2.4)



########################### Looking at optimal model combination

part1.12.optimal = readRDS("~/Dropbox/current_fluantigen/results/models/part1.twelvemonths.rds")

dataScaled.zero %>%
  mutate(above.zero = ifelse(three.months > 0, 1, 0)) -> dataScaled.zero

test.probs <- predict(part1.3.model, newdata=dataScaled.zero, type="response", allow.new.levels=TRUE)

test.ROC <- roc(predictor=test.probs, response=dataScaled.zero$above.zero)
test.ROC$auc

potential.thres = seq(.01, .99, .01)
compare.prediction = data.frame(test.probs, truth = as.factor(dataScaled.zero$above.zero), success = dataScaled.zero$success, current = dataScaled.zero$frequency)
predictedOutcomes = map(potential.thres, calculate_practicalResults, compare.prediction)
performance.statistics = map(predictedOutcomes, calculate_outcomeTypes, "class")


compare.prediction %>%
  mutate(above.thres = ifelse(test.probs > common.$thres, 1, 0)) %>%
  group_by(success, truth) %>%
  count(above.thres) 


optimal.three = data.frame(do.call("rbind", performance.statistics))
optimal.three$thres = potential.thres
optimal.three %>%
  arrange(desc(F1)) %>%
  slice(1) %>%
  select(thres) -> optimal.three.thres 


############# Comparing Optimal vs. Common Part 1 
common.three$model = "common";optimal.three$model = "optimal"
three.part1 = rbind(common.three, optimal.three); three.part1$month = "three"
common.six$model = "common"; optimal.six$model = "optimal"
six.part1 = rbind(common.six, optimal.six); six.part1$month = "six"
common.nine$model = "common"; optimal.nine$model = "optimal"
nine.part1 = rbind(common.nine, optimal.nine); nine.part1$month = "nine"
common.twelve$model = "common"; optimal.twelve$model = "optimal"
twelve.part1 = rbind(common.twelve, optimal.twelve); twelve.part1$month = "twelve"

part1.full.compare = rbind(three.part1, six.part1, nine.part1, twelve.part1)
part1.full.compare$month=factor(part1.full.compare$month, levels=c("three", "six", "nine", "twelve"),
                                labels=c("Three", "Six", "Nine", "Twelve"))

part1.full.compare %>%
  ggplot(aes(x = Sensitivity, y = Pos.Pred.Value, color = thres, shape = model)) +
  geom_point(size = 3) + scale_color_viridis() + facet_wrap(~month, nrow = 1) +
  labs(color = "Probability Threshold \n Cut-Off", shape = "Model Features") -> part1.plot

######################### Part 2
output.part2.3$month = "three"; output.part2.3$model = "common"
output.part2.6$month = "six"; output.part2.6$model = "common"
output.part2.9$month = "nine"; output.part2.9$model = "common"
output.part2.12$month = "twelve"; output.part2.12$model = "common"

part2.common = rbind(output.part2.3, output.part2.6, output.part2.9, output.part2.12)

output.part2.12.optimal$month = "twelve"; output.part2.12.optimal$model = "optimal"
output.part2.9.optimal$month = "nine"; output.part2.9.optimal$model = "optimal"
output.part2.6.optimal$month = "six"; output.part2.6.optimal$model = "optimal"
output.part2.3.optimal$month = "three"; output.part2.3.optimal$model = "optimal"

part2.optimal = rbind(output.part2.12.optimal, output.part2.9.optimal, output.part2.6.optimal, output.part2.3.optimal)

part2.full = rbind(part2.common, part2.optimal)
part2.full$month=factor(part2.full$month, levels = c("three", "six", "nine", "twelve"),
                        labels = c("Three","Six", "Nine", "Twelve"))
part2.full %>%
  ggplot(aes(x = pred^3, y = actual^3, color = model)) +
  geom_point(size = 3, alpha = .8) + facet_wrap(~month, nrow = 1) + 
  scale_color_viridis(discrete = TRUE) + labs(x = "Predicted Frequency", y = "Actual Frequency") + 
  geom_abline(intercept = 0, slope = 1) -> part2.plot.model

part2.full %>%
  ggplot(aes(x = pred^3, y = actual^3, color = success)) +
  geom_point(size = 3, alpha = .8) + facet_wrap(~month, nrow = 1) + 
  scale_color_viridis(discrete = TRUE, option = "B") + labs(x = "Predicted Frequency", y = "Actual Frequency") + 
  geom_abline(intercept = 0, slope = 1) 

write.csv(part2.full, file = "../results/models/frequency.predictions.csv")

combined.frequency.plot = plot_grid(part1.plot, part2.plot.model, nrow = 2)
save_plot(combined.frequency.plot, filename = "exploratory.figs/combined.frequency.plot.pdf",
          base_height = 8, base_aspect_ratio = 2.5)

########### Figure 4
key = c("0.01_Three", "0.56_Six", "0.54_Nine", "0.33_Twelve")

part1.full.compare %>%
  filter(model == "optimal") %>%
  mutate(data.key = paste0(thres, "_", month)) %>%
  filter(data.key %in% key) %>%
  ggplot(aes(x = Recall, y = Precision,color = month, label = month)) + geom_point(size = 5, alpha = .8) +
  scale_color_viridis(discrete = TRUE, option = "D") + 
  scale_x_reverse(limits = c(1, .95)) + 
  scale_y_continuous(limits = c(.83, 1.0)) + 
  geom_text(color = "black", nudge_x = .01, size = 8) + 
  labs(x = "Sensitivity", y = "Positive Predictive Value") +
  guides(color = FALSE) + 
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20)) -> part1.summary
part2.full %>%
  filter(model == "optimal") %>%
  ggplot(aes(x = pred^3, y = actual^3, color = success)) +
  geom_point(size = 3, alpha = .8) + facet_wrap(~month, nrow = 2) + 
  scale_color_manual(values = c("black", "grey")) + labs(x = "Predicted Frequency", y = "Actual Frequency", color = "Evolutionary Fate") + 
  geom_abline(intercept = 0, slope = 1) +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 18)) +
  theme(legend.text = element_text(size = 18),
        legend.title = element_text(size = 20)) +
  guides(color = guide_legend(keyheight = 2)) -> part2.summary


figure4 = plot_grid(part1.summary, part2.summary, rel_widths = c(1, 1.5))
save_plot(figure4, filename = "exploratory.figs/potential.figure4.pdf", base_height = 8, base_aspect_ratio = 2.1)
##################### Code Graveyard
#dependent.variable = "nine.months"
#folds.length = 5;
#dataPurged %<>% rename(name = trial) 
#folds = cvFolds(n = length(unique(dataPurged$name)), K = folds.length, type = "random")
#lmerRsquared=data.frame(matrix(nrow = folds.length, ncol = 2))
#lmerRmse =rep(0, folds.length)
#colnames(lmerRsquared) = c("marginal", "conditional")
#data.summary.mat = data.frame(matrix(nrow = folds.length, ncol  = 2))
#colnames(data.summary.mat) = c( "AIC","p.value")

#Initializing Variables 
#build = TRUE
#variable.list = c()
#variable.aic = c()
#variable.rmarginal = c()
#variable.rconditional = c()
#variable.rmse = c()
#anova.results = rep(0, folds.length)

#fixed.part.1 =  paste(dependent.variable, "~ value  + ")
#fixed.part.2 = " + (1 | name)"

#while(build == TRUE) {
# Step 1: Sets the formula and manipulates the data set
#if(length(variable.list) == 0 )  { # This is the first round
#  formula = paste(dependent.variable, "~value + (1 | name)") 
#  dataPurged %>%
#    ungroup() %>%
#    gather(key = variable, value = value, -one_of(dependent.variable), -name) %>%
#    mutate_at("value", as.numeric) -> data.scaled.l
#} else { # There are already significant variables 
#  variable.part = paste(variable.list, collapse = "+")
#  formula = paste0(fixed.part.1, variable.part, fixed.part.2)
#  formula.null = paste0(dependent.variable, "~", paste0(variable.list, collapse = "+"), fixed.part.2)
#  dataPurged %>%
#    ungroup() %>%
#    gather(key = variable, value = value, -one_of(dependent.variable), -name,  -one_of(variable.list)) %>%
#    mutate_at("value", as.numeric) -> data.scaled.l
#}
#model.results = ddply(data.scaled.l, .var = c("variable"), .fun = function(x) {
  # K-fold validation for each model based on a single term 
#  for(n in 1:folds.length) { 
#    test.ids = folds$subsets[folds$which==n]
#    test.trials = unique(data.scaled.l$name)[test.ids]
#    testdata = x[which(x$name %in% test.trials), ]
#    traindata = x[-which(x$name %in% test.trials),]
#    LMER <- lmer_alt(formula, data = traindata, return = "merMod")
#    data.summary.mat[n,] = c(glance(LMER)[3], tidy(LMER)[2,2])
#    lmerRsquared[n,] = r.squaredGLMM(LMER)
    # PRediction
#    lmer.probs <- predict(LMER, newdata=testdata, type="response", allow.new.levels=TRUE)
#    lmerRmse[n] <- rmse(actual = testdata[,dependent.variable], predicted = lmer.probs) 
#    if (length(variable.list) > 0) { # Compare new model against the null 
#      model.null = lme4::lmer(formula.null, data = traindata)
#      anova.test = anova(model.null, LMER, refit=FALSE)
#      anova.results[n] = anova.test$`Pr(>Chisq)`[2]
#    }
#  }
#  anova.mean = mean(anova.results)
#  rsquaredM=mean(lmerRsquared$marginal)
#  min.rsM = min(lmerRsquared$marginal)
#  max.rsM = max(lmerRsquared$marginal)
#  rsquaredC=mean(lmerRsquared$conditional)
#  min.rsC = min(lmerRsquared$conditional)
#  max.rscC = max(lmerRsquared$conditional)
#  rmse.mean = mean(lmerRmse)
#  data.summary = colMeans(as.data.frame(data.summary.mat))
#  return(t(c(rsquaredM = rsquaredM, min.rsM =min.rsM,max.rsM = max.rsM, 
#             rsquaredC = rsquaredC,min.rsC=min.rsC,max.rscC=max.rscC,
#             data.summary, anova = anova.mean, rmse.mean = rmse.mean)))# put the se of the performance 
#})
#variable = model.results %>% arrange(AIC) %>% slice(1)  # Select the variable that is the best
#if(length(variable.list) == 0 ) { # Is this the first? If yes, add this variable
#  variable.list = c(variable.list, variable$variable)
#  variable.aic = c(variable.aic, variable$AIC)
#  variable.rmarginal = c(variable.rmarginal, variable$rsquaredM)
#  variable.rconditional = c(variable.rconditional, variable$rsquaredC)
#  variable.rmse = c(variable.rmse, variable$rmse.mean)
#  print(paste0("Adding : ", variable$variable))
#} else if (variable$anova < .05) { # Determine if this new variable beats the null
#  variable.list = c(variable.list, variable$variable)
#  variable.rmarginal = c(variable.rmarginal, variable$rsquaredM)
#  variable.rconditional = c(variable.rconditional, variable$rsquaredC)
#  variable.aic = c(variable.aic, variable$AIC)
#  variable.rmse = c(variable.rmse, variable$rmse.mean)
#  print(paste0("Adding : ", variable$variable))
#} else { 
#  print(paste0("Final Single Term Model is: ", formula.null))
#  build = FALSE
#}
#}
#residual.test = lmer(formula.null, data = dataPurged)
#single.term.results = bind_cols(variable = variable.list, variable.rmarginal =  variable.rmarginal, 
#                                variable.rconditional = variable.rconditional, variable.aic = variable.aic,  variable.rmse = variable.rmse)
#write.csv(single.term.results, '../results/022818.sixmonths.arcsin.csv', row.names = FALSE)

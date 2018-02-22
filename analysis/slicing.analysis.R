##### Analysis Code for Random Point Surveillance 
rm(list=ls())
source('loadLibraries.R')
source('analysis_functions.R')
source('plotting_functions.R')


gather_data_emerge <- function(trial) {
  browser()
  trial.name = trial$name[1]
  console.file <- read.table(paste0(tropics.folder, trial.name, "/out.console.txt"), header = TRUE, fill = TRUE) %>% clean_console() 
  console.file %>%
    filter(postAntigen %in% trial$antigentype) %>%
    filter(!duplicated(postAntigen)) %>%
    select(-oriAntigen) %>%
    mutate_at("postAntigen", as.character) %>%
    right_join(trial, by = c("postAntigen" = "antigentype")) %>%
    mutate_at("day", as.character) -> trial
  
  track.antigen <- read.table(paste0(tropics.folder, trial.name, "/out.trackAntigenSeries.txt"), header = TRUE) %>%
    mutate_at(.vars = "day", as.character) 
  
  track.antigen %>%
    filter(day %in% trial$day) %>%
    right_join(trial, by = c("day" = "day")) -> trial
  
  viral.fitness.emergence <- read.table(paste0(tropics.folder,trial.name, "/out.viralFitnessSeries.txt"), header = TRUE) %>% 
    mutate_at(.var="day", as.character) 
  viral.fitness.emergence %>%
    filter(day %in% trial$day) %>%
    filter(!duplicated(day)) %>%
    select(-simDay) %>%
    right_join(trial, by = c("day" = "day")) -> trial
  
  antigen.frequencies[[eval(trial.name)]] %>%
    mutate_at("day", as.character) %>%
    filter(day %in% trial$day) %>%
    group_by(day) %>%
    slice(which.max(frequency)) %>%
    select(-antigentype) %>%
    right_join(trial, by  = c("day"="day")) -> trial
  
  ### function here for checking which ones aren't 
  return(trial)
}
gather_data_freq <- function(trial, surveillance.freq) {
  browser()
  trial.name = trial$name[1]
  surv.data = trial %>% select(postAntigen, success, name)
  
  # read in and extract information from antigen.frequencies for figuring out the first day 
  antigen.frequencies[[eval(trial.name)]] %>%
    filter(antigentype %in% trial$postAntigen) %>%
    group_by(antigentype) %>%
    filter(frequency > surveillance.freq) %>% 
    summarize(first.day = min(day)) %>%
    mutate_at("antigentype", as.character) %>%
    right_join(surv.data, by = c("antigentype" = "postAntigen")) -> surv.data
  # read in and extract from the viral fitness
  
  fitness <- read.table(paste0(tropics.folder, trial.name, "/out.viralFitnessSeries.txt"), header = TRUE)
  fitness %>%
    filter(day %in% surv.data$first.day) %>%
    filter(!duplicated(day)) %>%
    right_join(surv.data, by = c("day" = "first.day")) -> surv.data  ### at this point first.day becomes day 
  
  # read in and extract from the timeseries 
  timeseries = read.table(paste0(tropics.folder, trial.name, "/out.timeseries.txt"), header = TRUE)
  surv.data = adply(.data = surv.data, .margins = 1, function(surv) {
    timeseries %>%
      mutate(day.difference = abs(round(date*365)-surv$day)) %>%
      filter(day.difference <= 5) %>%
      slice(which.min(day.difference)) %>%
      select(diversity, tmrca, netau, serialInterval, antigenicDiversity, totalS, totalI, totalCases) %>%
      bind_cols(surv) -> surv
    return(surv)
  })
  
  surv.data = find_dominant_types_at_emerge(antigen.frequencies[[eval(trial.name)]]) %>%
    filter(day %in% surv.data$day) %>%
    rename(dominant.type = antigentype,  dominant.freq = frequency) %>%
    ungroup() %>%
    right_join(surv.data, by = "day")
  
  track.antigen <- read.table(paste0(tropics.folder, trial.name, "/out.trackAntigenSeries.txt"), header = TRUE)
  track.antigen %>%
    filter(day %in% surv.data$day) %>%
    select(-diversity, -tmrca, -netau, -serialInterval, -antigenicDiversity) %>%
    right_join(surv.data, by = "day") -> surv.data
  
  missing.trackAntigen = surv.data %>%
    filter(is.na(meanLoad))
  
  if(nrow(missing.trackAntigen) > 0) {
    surv.data %>%
      filter(!is.na(meanLoad)) -> surv.data
    missing.filled = adply(.data = missing.trackAntigen, .margins = 1, function(row) {
      track.antigen %>%
        dplyr::mutate(day.difference = abs(day-row$day)) %>%
        slice(which.min(day.difference)) %>%
        select(-diversity, -tmrca, -netau,-serialInterval,-antigenicDiversity, -day.difference) -> replacement.row
      return(replacement.row)
    })
    
    surv.data = rbind(surv.data, missing.filled)
  }
  return(surv.data)
}
gather_data_freq2 <- function(trial, surveillance.freq) {
  trial.name = trial$name[1]
  surv.data = trial %>% select(antigentype, success, name)
  
  # read in and extract information from antigen.frequencies for figuring out the first day 
  antigen.frequencies[[eval(trial.name)]] %>%
    filter(antigentype %in% trial$antigentype) %>%
    group_by(antigentype) %>%
    filter(frequency > surveillance.freq) %>% 
    summarize(first.day = min(day)) %>%
    mutate_at("antigentype", as.character) %>%
    right_join(surv.data, by = c("antigentype" = "antigentype")) -> surv.data
  # read in and extract from the viral fitness
  
  entropy[[eval(trial.name)]] %>%
    filter(day %in% surv.data$first.day) %>%
    filter(!duplicated(day)) %>%
    right_join(surv.data, by = c("day" = "first.day")) -> surv.data
  
  fitness <- read.table(paste0(tropics.folder, trial.name, "/out.viralFitnessSeries.txt"), header = TRUE)
  fitness %>%
    filter(day %in% surv.data$day) %>%
    filter(!duplicated(day)) %>%
    right_join(surv.data, by = "day") -> surv.data  ### at this point first.day becomes day 
  
  # read in and extract from the timeseries 
  timeseries = read.table(paste0(tropics.folder, trial.name, "/out.timeseries.txt"), header = TRUE)
  
  surv.data = adply(.data = surv.data, .margins = 1, function(surv) {
    timeseries %>%
      mutate(day.difference = abs(date - (surv$day/365-5))) %>%
      slice(which.min(day.difference)) %>%
      select(diversity, tmrca, netau, serialInterval, antigenicDiversity, totalS, totalI, totalCases) %>%
      bind_cols(surv) -> surv
    return(surv)
  })
  
  surv.data = find_dominant_types_at_emerge(antigen.frequencies[[eval(trial.name)]]) %>%
    filter(day %in% surv.data$day) %>%
    rename(dominant.type = antigentype,  dominant.freq = frequency) %>%
    ungroup() %>%
    right_join(surv.data, by = "day")
  
  #### There's a discrepancy here between which totalI I'm pulling in --- 
  track.antigen <- read.table(paste0(tropics.folder, trial.name, "/out.trackAntigenSeries.txt"), header = TRUE)
  track.antigen %>%
    filter(day %in% surv.data$day) %>%
    select(-diversity, -tmrca, -netau, -serialInterval, -antigenicDiversity) %>%
    right_join(surv.data, by = "day") -> surv.data
  
  
  viralTypeFitness <- read.table(paste0(tropics.folder, trial.name, "/out.typeViralFitness.txt"), header = TRUE)
  colnames(viralTypeFitness) = paste0("individual.", colnames(viralTypeFitness))
  surv.data = adply(.data = surv.data, .margins = 1, function(viralType) {
    viralTypeFitness %>%
      filter(individual.antigenType == viralType$antigentype) %>%
      filter(individual.day == viralType$day) %>%
      select(-individual.day, -individual.simDay, -individual.antigenType) %>%
      bind_cols(viralType) -> viralType
    return(viralType)
  })
  
  return(surv.data)
}
calculate_days_above <- function(antigen.frequencies, threshold) { 
  antigen.frequencies %>%
    group_by(antigentype) %>%
    dplyr::filter(frequency > threshold) %>%
    summarize(day.emerge = day[1],
              last.day = tail(day)[1],
              frequency.emerge = frequency[1],
              frequency.down = tail(frequency)[1]) %>%
    mutate(days = last.day-day.emerge) %>%
    dplyr::select(antigentype, days, day.emerge) %>%
    mutate_at("antigentype", as.character) 
}
replace_na_zeros = function(x) {
  #browser()
  already_lost = which(is.na(x$days))
  x$days[already_lost] = 0
  return(x)
}
determine_success_labels = function(x, max.rf) {
  x %>%
    mutate(success = ifelse(days > 45, "Est.", 
                            ifelse(final.max > max.rf, "Transient", "no"))) %>%
    mutate(day.success = ifelse(success == "Est.", day.emerge + 45,  NA))
}
filter_out_loss = function(x) {
  x %>%
    filter(success != "no")
}
clean_gathered_data = function(entry) {
  trial.name = entry$name
  
  track.antigen <- read.table(paste0(tropics.folder, trial.name, "/out.trackAntigenSeries.txt"), header = TRUE) 
  entry.day = as.integer(entry$day)
  
  #Determine alternate day 
  track.antigen %>%
    mutate(day.off = abs(day-entry.day)) %>%
    slice(which.min(day.off)) %>%
    select(-day.off) -> selected.day 
  entry$day = selected.day$day
  
  ## Get rid of ones that are NA
  na.indices = sapply(entry, is.na)
  entry %>%
    select(-which(na.indices == TRUE)) %>%
    left_join(selected.day) %>%
    mutate_at("day", as.character) -> entry
  
  if(ncol(entry) == 32) {
    return(entry)
  } else {
    viral.fitness.emergence <- read.table(paste0(tropics.folder,trial.name, "/out.viralFitnessSeries.txt"), header = TRUE) %>% 
      mutate_at(.var="day", as.character) 
    viral.fitness.emergence %>%
      filter(day %in% entry$day) %>%
      filter(!duplicated(day)) %>%
      select(-simDay) %>%
      left_join(entry, by = c("day" = "day")) -> entry
    
    antigen.frequencies[[eval(trial.name)]] %>%
      mutate_at("day", as.character) %>%
      filter(day %in% entry$day) %>%
      group_by(day) %>%
      slice(which.max(frequency)) %>%
      select(-antigentype) %>%
      right_join(entry, by  = c("day"="day")) -> entry
    return(entry)
  }
}
calculate_ratios = function(data.set) {  
  data.set %>%
    mutate(ratio.mutation = individual.meanMut/meanLoad,
           ratio.meanR = exp(individual.meanR)/meanR,
           ratio.varR = individual.varR/varR,
           ratio.meanBeta = exp(individual.meanBeta)/meanBeta,
           ratio.varBeta = individual.varBeta/varBeta,
           ratio.meanSigma = individual.meanSigma/meanSigma,
           ratio.varSigma = individual.varSigma/varSigma) 
}
remove_columns = function(data.set) {
  excluded.variables = c("N", "S", "I", "R", "cases", "cumulativeTypes", "dominant.type", "totalCases")
  data.set %>%
    select(-one_of(excluded.variables)) 
}
calculate_entropy = function(antigen.frequency) {
  antigen.frequency %>%
    group_by(day) %>%
    mutate(log.frequency = log(1/frequency),
           ind.term = frequency * log.frequency) %>%
    summarize(entropy = sum(ind.term))
}


#### new functions
collect_poplevel = function(time, trial.name) {
  
  fitness <- read.table(paste0(tropics.folder, trial.name, "/out.viralFitnessSeries.txt"), header = TRUE)
  fitness.population = fitness %>% filter(day == time)
  date = as.data.frame(fitness$day)
  
  timeseries = read.table(paste0(tropics.folder, trial.name, "/out.timeseries.txt"), header = TRUE)
  timeseries.population = timeseries %>% bind_cols(date) %>%
    filter(fitness$day == time) %>%
    rename(day = 'fitness$day') %>%
    select(date, diversity, tmrca, netau, serialInterval, antigenicDiversity, totalS, totalI, day)
  
  track.antigen <- read.table(paste0(tropics.folder, trial.name, "/out.trackAntigenSeries.txt"), header = TRUE)
  antigen.population = track.antigen %>% filter(day == time) %>% select(day, meanLoad, antigenicTypes)
  
  entropy.population = entropy[[eval(trial.name)]] %>% filter(day == time)
  
  pop.level.list = list(fitness.population, timeseries.population, antigen.population, entropy.population)
  join_all(pop.level.list, by = "day", type = "left")
} #### THIS IS GOOD 
find_antigens_present <- function(time, trial.name) {
  trial.antigen.frequencies = antigen.frequencies[[eval(trial.name)]]
  trial.antigen.frequencies %>%
    filter(day == time) %>%
    select(day, antigentype, frequency)
}
create_trial_dataset <- function(trial.name, sample.times) { 
  pop.level.list = map(sample.times, collect_poplevel, trial.name)
  pop.level.df = do.call("rbind", pop.level.list) %>% arrange(day)
  
  
  ### Identify which antigens are present at these times 
  antigen.present.list = map(sample.times, find_antigens_present, trial.name)
  antigen.present.df = do.call("rbind", antigen.present.list) %>% arrange(day)
  
  
  ### Filter out those that are above the threshold of capture (1%) have already been successful
  antigen.present.df %>%
    filter(frequency > .01) -> antigen.present.captured
  
  ### Filter out ones that have reached their peak already 
  subset.analyze %>%
    filter(name == trial.name) %>%
    filter(antigentype %in% antigen.present.captured$antigentype) %>%
    select(antigentype, day.success) %>%
    mutate_at("antigentype", as.factor) -> antigen.subset.emerge
  
  antigen.present.captured %>%
    mutate_at(.vars = "antigentype", as.factor) %>%
    arrange(antigentype) %>%
    left_join(antigen.subset.emerge, by = "antigentype") %>%
    mutate(stage = ifelse(day.success < day, "peaked", "rising")) %>%
    filter(stage != "peaked" | is.na(stage)) %>%
    mutate(stage  = "rising") -> antigen.present.analyze
  
  ### at this point ahve to go in and get their info
  
  viralTypeFitness <- read.table(paste0(tropics.folder, trial.name, "/out.typeViralFitness.txt"), header = TRUE)
  colnames(viralTypeFitness) = paste0("individual.", colnames(viralTypeFitness))
  
  individual.fitness = ddply(antigen.present.analyze,.variables = "day", function(x) {
    viralTypeFitness %>%
      filter(individual.day == x$day[1]) -> viralType.day
    viralType.day %>%
      filter(individual.antigenType %in% x$antigentype) %>%
      select(-individual.simDay) %>%
      rename(day = individual.day) %>%
      mutate_at(.vars = "individual.antigenType", as.character)
  })  
  
  individual.fitness = left_join(individual.fitness,antigen.present.analyze, by = c("day" = "day", "individual.antigenType" = "antigentype"))
  ############### Combine the population and individual data
  combined.data = left_join(individual.fitness, pop.level.df) 
  combined.data = combined.data %>% calculate_ratios()
  
  ############# combine here ultimate success or not 
  full.data[[eval(trial.name)]]  %>%
    filter(antigentype %in% combined.data$individual.antigenType) %>%
    select(antigentype, success) %>%
    right_join(combined.data, by = c("antigentype"="individual.antigenType")) %>%
    filter(antigentype != "0") -> combined.data
  
  combined.data %>%
    gather(key = variable, value = value, -day, -antigentype, -success) %>%
    mutate(trial = trial.name) -> combined.data.l
  
  return(combined.data.l)
}
generate_second_time = function(sample.time) {
  possible.second = seq(from = sample.time+14, to = sample.time + 4*14, by = 14)
  second.time = sample(possible.second, 1)
  return(second.time = second.time)
} ######### new function 
apply_group_number <- function(data) {
  data %>%
    mutate(., group = group_indices(.,day))
}
create_time_variable <- function(data) { 
  data %>%
    apply_group_number() %>%
    select(antigentype, success, day, trial, group) %>%
    distinct(antigentype, success, day, trial, group) %>%
    gather(key = variable, value = value, day) -> time.variable
  data %>%
    apply_group_number() %>%
    select(-day) %>%
    rbind(time.variable)
}
calculate_diff <- function(data) { 
  data %>%
    filter(!(variable %in% remove.variables)) %>%
    mutate_at("value", as.numeric) %>% 
    spread(key = time.stamp, value = value, fill = 0) %>%
    arrange(antigentype, group) %>% 
    mutate(diff = t2-t1) %>%
    gather(key = time, value = value, t1, t2, diff) 
}


###########

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
rel.frequency.thres = .01 # Threshold at which we start monitoring/exepct to pick up in surveillance 
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
possible.sample = seq(from = 1833, to = 10947, by = 14)
sample.time.one = sample(possible.sample, 50, replace=FALSE) 
sample.time.two = adply(sample.time.one, .margins = 1, generate_second_time,.id = NULL)
sample.times = data.frame(t1 = sample.time.one, t2= sample.time.two$V1) %>% arrange(t1)
head(sample.times)

######################## Function to collect population level at those times
 # For trial time 1 
data.1.list = map(trial.dirs, create_trial_dataset, sample.times$t1)
data.1.df = map(data.1.list, create_time_variable) %>% 
  map(function(x) x %>% mutate(time.stamp = 't1')) %>%
  dplyr::bind_rows(.id = 'list')

data.1.df %>%
  distinct(antigentype,success,trial,group) %>%
  group_by(success, group) %>%
  summarize(n = n()) %>%
  group_by(group) %>%
  mutate(rel.freq  = n/sum(n))%>%
  arrange(group)

remove.variables = c("date", "day", "day.success", "stage", "simDay")
data.1.df %<>%
  filter(!(variable %in% remove.variables))

###############################
#data.2.list = map(trial.dirs, create_trial_dataset, sample.times$t2)
#data.2.df = map(data.2.list, create_time_variable) %>% 
#  map(function(x) x %>% mutate(time.stamp = 't2')) %>%
#  dplyr::bind_rows(.id = 'list')

#data.time.combined = rbind(data.1.df, data.2.df)

################## Need to this for every trial 
# Back into list format based on trial 
#full.data.list = dlply(data.time.combined, "trial", identity)

#remove.variables = c("day.success", "stage")

## Calculate diff for each trial   
#list.calculate.diff = map(full.data.list, calculate_diff)

## Combine 
#full.data.df = do.call("rbind", list.calculate.diff)
#full.data.df %>% 
#  filter(variable == "day" & success == "Est.") %>% arrange(antigentype, group) %>%
#  filter(time == "t2" & value == 0)

# Calculate diff in sample.times
#sample.times %<>% mutate(time.diff = t2-t1, group = row_number())

#sample.times %>%
#  select(time.diff, group) %>%
#  left_join(full.data.df) %>%
#  head()

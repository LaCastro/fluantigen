#### If only looking at transient -- should identify which ones those are right away before 
#### collecting the rest of the data on them
rm(list=ls())
source('loadLibraries.R')
source('analysis_functions.R')
source('plotting_functions.R')

gather_data_emerge <- function(trial) {
  
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
      mutate(day.difference = abs(round(date*365)-surv$day)) %>%
     # filter(date == (surv$day/365)) %>%
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
    dplyr::select(antigentype, days) %>%
    mutate_at("antigentype", as.character) 
}
replace_na_zeros = function(x) {
  #browser()
  already_lost = which(is.na(x$days))
  x$days[already_lost] = 0
  return(x)
}
determine_success_labels = function(x) {
  x %>%
    mutate(success = ifelse(days > 45, "Est.", 
                            ifelse(final.max > .02, "Transient", "no")))
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


###########
# Step 1. Put together list trials of transient and successful antigens ; this will be based on max.freq, and days above
tropics.folder = "../data/tropics_30/eligible/"
tropics.folder = "../data/tropics/eligible/"
trial.dirs = dir(tropics.folder)

# This is a slow step
antigen.frequencies = map(trial.dirs, function(x) read.table(paste0(tropics.folder,x,"/out.antigenFrequencies.txt"), header = TRUE)) 
names(antigen.frequencies) = trial.dirs

######## Step 2: For each list entry need to calculate max frequency and days above
thres = .2

# Calculating the two criteria for determining success 
antigen.frequencies %>%
  map(find_max_frequency) %>%
  map(function(x) mutate_at(x, "antigentype", as.character)) -> max.frequencies

antigen.frequencies %>%
  map(function(x) calculate_days_above(x,threshold = .2)) -> days.above

# Creating the combined data frame and using the criteria to assign a label 
full.data = map2(max.frequencies, days.above, left_join) %>%
  map(replace_na_zeros) %>% 
  map(determine_success_labels)
#rm(max.frequencies, days.above)

full.data.df = do.call("rbind", full.data)
full.data.analyze = full.data.df %>% filter(antigentype != 0)

full.data.analyze %>%
  group_by(success) %>%
  summarize(num = n()) %>%
  mutate(freq = num/sum(num))

####### Step 3: Subsetting for just transient and yes; and getting rid of zeros 
subset.data = map(full.data, filter_out_loss) 
subset.df = do.call("rbind", subset.data)
subset.df$name = rep(trial.dirs, sapply(subset.data, nrow))

subset.df %>%
  filter(antigentype != 0) -> subset.analyze

subset.analyze %>%
  group_by(success) %>%
  


####### Step 4: Create Meta Data for these 
#meta.data = ddply(.data = subset.analyze, .variables = "name", function(trial) gather_data_emerge(trial)) %>%
#  mutate(id = paste0(postAntigen, "_", name))
#meta.data %>%
#  filter(is.na(meanLoad) | is.na(frequency) | is.na(varBeta)) -> missing
#missing.filled = adply(.data = missing, .margins = 1, clean_gathered_data)
#meta.data %>%
#  filter(!(id %in% missing$id)) %>%
#  bind_rows(missing.filled) -> meta.data


######### Step 4 - Calculate Entropy
entropy = map(antigen.frequencies, calculate_entropy)
names(entropy) = trial.dirs




####### Step 5: Create Meta Data at different time points
freq.one = ddply(.data = subset.analyze, .variables = "name", function(trial) gather_data_freq2(trial, surveillance.freq = .01))
freq.two = ddply(.data = subset.analyze, .variables = "name", function(trial) gather_data_freq2(trial,surveillance.freq = .02))
freq.three = ddply(.data = subset.analyze, .variables = "name", function(trial) gather_data_freq2(trial,surveillance.freq = .03))
freq.four = ddply(.data = subset.analyze, .variables = "name", function(trial) gather_data_freq2(trial,surveillance.freq = .04))
freq.five = ddply(.data = subset.analyze, .variables = "name", function(trial) gather_data_freq2(trial, surveillance.freq = .05))

freq.list = list(freq.one, freq.two, freq.three, freq.four, freq.five) 

######## Step 6: Remove extra rows 
freq.list = map(freq.list, remove_columns)

######## Step 7: Calculate Ratio of fitness 
freq.list.ratio = map(freq.list, calculate_ratios)

###### Step 8: Calculate Difference Data Sets 
freq.df =  do.call("rbind", freq.list.ratio)
freq.df$freq = rep(c("freq.01", "freq.02", "freq.03", "freq.04", "freq.05"), sapply(freq.list.ratio, nrow))

#### Step 8.5 - remove those that span  the burn in period 
freq.df %>%
  group_by(antigentype, name) %>%
  filter(day == 1833) %>%
  summarize(antigen.to.exclude = unique(antigentype))%>%
  mutate(.id = paste0(name, "_", antigen.to.exclude)) -> not.full.dataset

freq.df %>%
  mutate(.id = paste0(name, "_", antigentype)) %>%
  filter(!(.id %in% not.full.dataset$.id)) %>%
  select(-.id) -> freq.df.subset

############## Section 9 - Calculate Growth Differences 

calculate_growth_phase = function(x, phase1, phase2) {
  (x[phase2]-x[phase1])
}

freq.df.subset %>%
  group_by(antigentype, name) %>%
  gather(key = variable, value = value, -freq,-antigentype,-success,-name) %>%
  arrange(antigentype) %>%
  group_by(name, antigentype, variable) %>%
  summarize(first.gp = calculate_growth_phase(value, phase1 = 1, phase2 = 3),
            second.gp = calculate_growth_phase(value, phase1 = 3, phase2 = 5),
            success = success[1]) %>%
  ungroup() %>%
  mutate_at("name", as.factor) -> diff.df

diff.df %>%
  select(-second.gp) %>%
  spread(key = variable, value = first.gp) -> growth.1.subset

################ Step 9: Calculate Acceleration
diff.df %>%
  mutate(accelerate = second.gp - first.gp) -> diff.df

diff.df %>%
  select(-first.gp, -second.gp) %>%
  spread(key = variable, value = accelerate) -> accelerate


 

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
                            ifelse(final.max > .1, "Transient", "no")))
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

# Step 1. Put together list trials of transient and successful antigens ; this will be based on max.freq, and days above
tropics.folder = "../data/tropics/tropics_20yr/"
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
names(full.data) = trial.dirs

####### Step 3: Subsetting for just transient and yes; and getting rid of zeros 
subset.data = map(full.data, filter_out_loss) 
subset.df = do.call("rbind", subset.data)
subset.df$name = rep(trial.dirs, sapply(subset.data, nrow))

subset.df %>%
  filter(antigentype != 0) -> subset.analyze

####### Step 4: Create Meta Data for these 
meta.data = ddply(.data = subset.analyze, .variables = "name", function(trial) gather_data_emerge(trial))
meta.data = mutate(meta.data, id = paste0(postAntigen, "_", name))
meta.data %>%
  filter(is.na(meanLoad) | is.na(frequency) | is.na(varBeta)) -> missing
meta.data %>%
  filter(!(id %in% missing.filled$id)) -> meta.data
missing.filled = adply(.data = missing, .margins = 1, clean_gathered_data)
meta.data = rbind(meta.data, missing.filled)


####### Step 5: Create Meta Data at different time points
meta.data %>% filter(name == "tropics_100") -> trial

gather_data_freq(trial, surveillance.freq) {
  trial.name = trial$name[1]
  surv.data = trial %>% select(postAntigen, success, name)
  
  trial.meta.data %>%
    filter(success == eval(type)) %>%
    distinct(postAntigen) -> selected.antigens
  
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
    right_join(surv.data, by = c("day" = "first.day")) -> surv.data
  
  # read in and extract from the timeseries 
  timeseries = read.table(paste0(tropics.folder, trial.name, "/out.timeseries.txt"), header = TRUE)
    # Gonna need the timeseries modifications in here 
  surv.data = adply(.data = surv.data, .margins = 1, function(surv) {
     timeseries %>%
        mutate(day.difference = abs(round(date*365)-surv$day)) %>%
        filter(day.difference <= 5) %>%
        slice(which.min(day.difference)) %>%
        select(diversity, tmrca, netau, serialInterval, antigenicDiversity, totalS, totalI, totalCases) %>%
        bind_cols(surv) -> surv
    return(surv)
    })
    
############# GOT TO THIS POINT ! 
  
    ## find dominant type 
    dominant.types = find_dominant_types_at_emerge(antigen.frequencies) %>%
      filter(day %in% full.freq.summary$first.day) %>%
      rename(dominant.type = antigentype,  dominant.freq = frequency) %>%
      ungroup() %>%
      mutate_at(.vars = "day", as.character)
    
    full.freq.summary %>%
      left_join(dominant.types, by = c("first.day" = "day")) -> full.freq.summary
    
    # read in and extract track antigen closest to that day
    track.antigen <- read.table(paste0(sim.dir, "/out.trackAntigenSeries.txt"), header = TRUE)
    
    track.antigen.summary = adply(.data = day.freq.selected, .margins = 1, function(antigen) {
      track.antigen %>%
        mutate(day.difference = abs(day-antigen$first.day)) %>%
        filter(day.difference == min(day.difference)) %>%
        dplyr::select(1, 12:14) -> track.antigen.summary
      return(track.antigen.summary)
    })
    ### find the closest to the day 
    track.antigen.summary$first.day = as.character(track.antigen.summary$first.day)
    
    full.freq.summary %>%
      left_join(track.antigen.summary, by = c("first.day", "postAntigen")) %>%
      rename(-> full.freq.summary
             
             colnames(full.freq.summary)[2] = "day"
             return(full.freq.summary)
}



find_data_at_freq <- function(sim.dir, trial.meta.data, surveillance.freq, type) {
  # filter data to look at either successful or unsuccesful depending on desire
  trial.meta.data %>%
    filter(success == eval(type)) %>%
    distinct(postAntigen) -> selected.antigens
  
  # read in and extract information from antigen.frequencies 
  antigen.frequencies <- read.table(paste0(sim.dir, "/out.antigenFrequencies.txt"), header = TRUE)
  
  # First Day that antigen is above surveillance threshold 
  day.freq.selected = adply(.data = selected.antigens, .margins = 1, function(antigen) {
    antigen.frequencies %>%
      filter(antigentype == antigen$postAntigen) %>%
      filter(frequency > surveillance.freq) %>%
      summarize(first.day = min(day))
  }) 
  
  # If there are no antigens that reached threshold but were not successful, move on 
  if(nrow(day.freq.selected) ==0) {
    return()
  } else { 
    # read in and extract from the viral fitness
    fitness <- read.table(paste0(sim.dir, "/out.viralFitnessSeries.txt"), header = TRUE)
    fit.summary = adply(.data = day.freq.selected, .margins = 1, function(antigen) {
      fitness %>%
        filter(day == antigen$first.day) %>%
        dplyr::select(-day) -> state
      antigen.state = cbind(antigen, state)
    })
    
    # read in and extract from the timeseries 
    timeseries = read.table(paste0(sim.dir, "/out.timeseries.txt"), header = TRUE)
    
    # Gonna need the timeseries modifications in here 
    timeseries.summary = adply(.data = day.freq.selected, .margins = 1, function(antigen) {
      timeseries %>%
        mutate(day.difference = date*365-antigen$first.day) %>% 
        filter(abs(day.difference) <= 5) %>% # helps link up the date and day
        dplyr::select(1:11)
    })
    
    left_join(fit.summary, timeseries.summary, by = c("postAntigen", "first.day"))  -> full.freq.summary
    full.freq.summary$first.day = as.character(full.freq.summary$first.day)
    
    ## find dominant type 
    dominant.types = find_dominant_types_at_emerge(antigen.frequencies) %>%
      filter(day %in% full.freq.summary$first.day) %>%
      rename(dominant.type = antigentype,  dominant.freq = frequency) %>%
      ungroup() %>%
      mutate_at(.vars = "day", as.character)
    
    full.freq.summary %>%
      left_join(dominant.types, by = c("first.day" = "day")) -> full.freq.summary
    
    # read in and extract track antigen closest to that day
    track.antigen <- read.table(paste0(sim.dir, "/out.trackAntigenSeries.txt"), header = TRUE)
    
    track.antigen.summary = adply(.data = day.freq.selected, .margins = 1, function(antigen) {
      track.antigen %>%
        mutate(day.difference = abs(day-antigen$first.day)) %>%
        filter(day.difference == min(day.difference)) %>%
        dplyr::select(1, 12:14) -> track.antigen.summary
      return(track.antigen.summary)
    })
    ### find the closest to the day 
    track.antigen.summary$first.day = as.character(track.antigen.summary$first.day)
    
    full.freq.summary %>%
      left_join(track.antigen.summary, by = c("first.day", "postAntigen")) %>%
      rename(-> full.freq.summary
             
             colnames(full.freq.summary)[2] = "day"
             return(full.freq.summary)
  }
}
find_data_at_freq_all <- function(dir, correct.trials, meta.data, surveillance.freq, type) {
  browser()
  meta.data %>%
    dplyr::select(-cases, -simDay, -dominant.type) %>%
    filter(final.max > surveillance.freq) %>%  
    gather(key = metric, value = emergence.value, -.id, -postAntigen, -success) %>%
    mutate(postAntigen = as.numeric(postAntigen)) %>%
    arrange(postAntigen, .id) %>%
    filter(metric != "N") -> tidy.data 
  
  freq.data <- lapply(correct.trials, function(trial) {
    tidy.data %>%
      filter(.id == trial) -> trial.tidy.data
    
    meta.data = find_data_at_freq(sim.dir = paste0(dir,trial), 
                                  surveillance.freq = surveillance.freq, 
                                  trial.meta.data = trial.tidy.data,
                                  type = type)
    meta.data$.id = trial
    return(meta.data)
  })
  
  empty.lists=which(sapply(freq.data,is.null))
  if(length(empty.lists > 0)) { 
    freq.data = freq.data[-empty.lists]
  }
  freq.data.all = rbindlist(freq.data, fill = TRUE)
  return(freq.data.all)
}






###Calculates the number of days an antigen type is present above a certain frequency


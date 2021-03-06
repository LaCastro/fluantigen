rm(list=ls())
source('analysis_functions.R')
source('plotting_functions.R')

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

data_at_freq <- function(dir, correct.trials, surveillance.freq, meta.data, summary.infection) {
  
  freq.no = find_data_at_freq_all(dir = dir, correct.trials = correct.trials,
                                  surveillance.freq =  surveillance.freq, meta.data = meta.data, type = "no")
  freq.yes = find_data_at_freq_all(dir = dir, correct.trials =correct.trials,
                                   surveillance.freq = surveillance.freq, meta.data = meta.data, type = "yes")
  
  freq.both = rbind(data.frame(success = "yes", freq.yes),
                    data.frame(success = "no", freq.no))
  
  
  left_join(freq.both, summary.infection) %>%
    mutate(normalize.I = (infected-min.I)/(max.I-min.I)) -> freq.both
  return(freq.both)
}

############################################################################################################
data.folder = "../data/tropics/tropics_20yr"

antigen.data = create_meta_data_all(dir = data.folder)
timeseries = read_outputfiles(data.folder, "/out.timeseries.txt")

### For north, going to calculate the timeseries relative to what it's like below 10 years
infected.range = calculate_max_infected(timeseries)
antigen.data = normalize_infection(meta.data = antigen.data, infected.range)

thres = .2
antigen.frequencies = read_outputfiles(data.folder, "/out.antigenFrequencies.txt")

## this could eventually be re-done so the calculate_days_above is happening when generating the data set 
days.above.thres = calculate_days_above_thres(antigen.frequencies, threshold = thres)
antigen.data %>% left_join(days.above.thres, by = c("postAntigen" = "antigentype", ".id" = ".id")) -> antigen.data

already_lost = which(is.na(antigen.data$days.above))
antigen.data$days.above[already_lost] = 0

antigen.data %>%
  mutate(success = ifelse(days.above > 45, "Est.", ifelse(final.max > .1, "Transient", "no"))) -> antigen.data

antigen.data %>%
  group_by(.id) %>%
  summarize(max.I = max(max.I),
            min.I = min(min.I)) -> infection.summary 

## have to make data tidy before this 
correct.trials = unique(antigen.data$.id)

freq.02 = data_at_freq(dir=data.folder, meta.data = antigen.data, correct.trials = correct.trials,
                         surveillance.freq = .02, summary.infection=infection.summary)
freq.15 = data_at_freq(dir=data.folder, meta.data = antigen.data, correct.trials = correct.trials,
                         surveillance.freq = .15, summary.infection=infection.summary)


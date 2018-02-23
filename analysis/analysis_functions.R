###### Functions

################### Data Cleaning and Data Prep ################
read_outputfiles <- function(dir, type ) {
  # Read in output file type from directory 
  # Reads in all trials in that directory 
  file.list = list.dirs(dir, full.names = FALSE, recursive = FALSE)
  population.data <- lapply(file.list, function(.file) {
    output.file <- read.table(paste0(dir,.file, type), header = TRUE)
    print(.file) # when testing if the formatting is correct 
    output.file
  })
  names(population.data) = file.list
  population.data = rbindlist(population.data, idcol = TRUE)
  return(population.data)
}

read_console_files <- function(dir) {
  # Reads console output files
  # Reads in all trials in that directory 
  file.list = list.dirs(dir, full.names = FALSE, recursive = FALSE)
  population.data <- lapply(file.list, function(.file) {
    output.file <- read.table(paste0(dir,.file, "/out.console.txt"), header = TRUE, fill = TRUE)
    # Clean console file here 
    output.file = clean_console(output.file)
    return(output.file)
  })
  names(population.data) = file.list
  #population.data = rbindlist(population.data, idcol = TRUE)
  return(population.data)
}

clean_console <- function(console.file) {
  # Cleans up console-gets rid of last two lines and if line break wasn't put in 
  output.file = console.file[1:(nrow(console.file)-2), ]
  na.rows = output.file[is.na(output.file$postAntigen),]
  na.indicies = which(is.na(output.file$postAntigen))
  modified.rows = output.file[na.indicies-1,]
  if (length(na.indicies)>0) {
    for(i in 1:length(na.indicies)) {
      index = na.indicies[i]
      output.file[index, 2:6] = unname(na.rows[i,1:5])
      output.file[index, 1] = modified.rows[i,6]
      distance =  na.rows[i, 3]
      output.file[index, 4] = as.character(distance)
      output.file[index, 1] = modified.rows[i,6]
    }
    output.file = output.file[-(na.indicies-1), ]
  }
  return(output.file)
}

remove_trials_data <- function(data, trials) {
  # Remove trials from meta data if dynamics 
  # do not reflect empirical 
  data %>%
    filter(.id %in% trials)
}

find_antigen_emergence <- function(antigen.data){
  # Find entries of novel antigens 
  antigen.data %>%
    filter(!duplicated(postAntigen)) %>%
    filter(postAntigen != 0) -> novel.types 
  return(novel.types)
}



find_dominant_types_at_emerge <- function(antigen.frequencies) {
  # Find the frequency of the dominant type circulating at emergence 
  antigen.frequencies %>%
    group_by(day) %>%
    slice(which.max(frequency)) 
}


find_max_frequency <- function(antigen.frequencies) {
  # max frequency antigen type achieved 
  antigen.frequencies %>%
    group_by(antigentype) %>%
    summarize(final.max = max(frequency))
}

create_meta_data <- function(sim.dir) {
  browser()
  # Combines all the output files for novel antigens
  # reading in and cleaning up the console file 
  console.file <- read.table(paste0(sim.dir, "/out.console.txt"), header = TRUE, fill = TRUE) %>% clean_console()

  # find the time when novel antigens emerge
  novel.types = find_antigen_emergence(console.file) %>% mutate_at("day", as.character)

  # read in output file that tracks state of simulation when antigens emerge
  track.antigen <- read.table(paste0(sim.dir, "/out.trackAntigenSeries.txt"), header = TRUE) %>%
    mutate_at(.vars = "day", as.character)
  meta.data  = left_join(x = novel.types, y = track.antigen, "day") 
  
  # record the exact days when new novel types are generated 
  days.of.emergence = meta.data$day
  
  ## Read in viral fitness, dplyr::select on the days that correspond to emergence 
  viral.fitness.emergence <- read.table(paste0(sim.dir, "/out.viralFitnessSeries.txt"), header = TRUE) %>%
    filter(day %in% days.of.emergence) %>%
    mutate(day = as.character(day)) %>%
    dplyr::select(-simDay) %>%
    filter(!duplicated(day)) 
 
   meta.data %>%
    left_join(viral.fitness.emergence, by = "day") -> meta.data
  
  ## combine dominant frequency circulating at time of emergence
  antigen.frequencies <- read.table(paste0(sim.dir, "/out.antigenFrequencies.txt"), header = TRUE) 
   dominant.types = find_dominant_types_at_emerge(antigen.frequencies) %>%
    filter(day %in% days.of.emergence) %>%
    dplyr::rename(dominant.type = antigentype, dominant.freq = frequency) %>%
    ungroup() %>%
    mutate_at(.vars = "day", as.character)
  
  meta.data %>%
    left_join(dominant.types, by = "day") -> meta.data
  
  # combine maximum frequency the strain itself ever achieved
  maximum.freq.type <- find_max_frequency(antigen.frequencies) %>% mutate_at("antigentype", as.character)
  meta.data %>%
    mutate_at("postAntigen", as.character) %>%
    left_join(maximum.freq.type, by = c("postAntigen" = "antigentype")) -> meta.data
  
  life.span <- calculate_total_life(antigen.frequencies) %>% mutate_at("antigentype", as.character)
  meta.data %>% left_join(life.span, by = c("postAntigen" = "antigentype")) -> meta.data
  
  # Differentiate whether it was sucessful or not
  #days.above.thres = _above_thres(antigen.frequencies, success.criteria$threshold)
  #days.above.thres$antigentype = as.character(days.above.thres$antigentype)
  #meta.data %>% left_join(days.above.threshold, by = c("postAntigen" = "antigentype")) -> meta.data
  #meta.data %>%
  #  mutate(success = ifelse((final.max > success.criteria$freq & life.length > success.criteria$length.days), "yes", "no")) -> meta.data
  
  return(meta.data)
}

create_meta_data_all <- function(dir) {
  # Runs create_meta_data on all trials in a directory 
  file.list = list.dirs(dir, full.names = FALSE, recursive = FALSE)
  population.data <- lapply(file.list, function(.file) {
    meta.data = create_meta_data(paste0(dir, sim.dir = .file))
    return(meta.data)
  })
  names(population.data) = file.list
  population.data = rbindlist(population.data, idcol = TRUE)
  return(population.data)
}



create_short_data <- function(sim.dir) {
  ### Combines all the output files for novel antigens
  # reading in and cleaning up the console file 
  
  console.file <- read.table(paste0(sim.dir, "/out.console.txt"), header = TRUE, fill = TRUE) %>% clean_console()
  # find the time when novel antigens emerge
  novel.types = find_antigen_emergence(console.file) %>% mutate_at("postAntigen", as.character)
  
  ## combine dominant frequency circulating at time of emergence
  antigen.frequencies <- read.table(paste0(sim.dir, "/out.antigenFrequencies.txt"), header = TRUE)
  life.span <- calculate_total_life(antigen.frequencies) %>% mutate_at("antigentype", as.character)
  
  # combine maximum frequency the strain itself ever achieved
  maximum.freq.type <- find_max_frequency(antigen.frequencies) %>%
    mutate_at("antigentype", as.character)
 
  novel.types %>% 
    left_join(life.span, by = c("postAntigen" = "antigentype")) %>%
    left_join(maximum.freq.type, by = c("postAntigen" = "antigentype")) -> novel.types
              
  return(novel.types)
}

create_short_data_all <- function(dir) {
  # Runs create_meta_data on all trials in a directory 
  file.list = list.dirs(dir, full.names = FALSE, recursive = FALSE)
  population.data <- lapply(file.list, function(.file) {
    meta.data = create_short_data(paste0(dir, sim.dir = .file))
    return(meta.data)
  })
  names(population.data) = file.list
  population.data = rbindlist(population.data, idcol = TRUE)
  return(population.data)
}

count_n_antigens <- function(antigen.record) {
  # Return number of unique antigens in a record for each trial
  antigen.record %>%
    group_by(.id) %>%
    summarize(unique.antigens = n_distinct(antigentype)) 
}

count_success_antigens <- function(antigen.record, threshold.freq, threshold.day) {
  # Counnt the number of antigens that meet a frequency and lifespan criteria 
  antigen.record %>%
    group_by(.id, antigentype) %>%
    filter(frequency > threshold.freq) %>%
    summarize(days.alive = day[n()]-day[1]) %>%
    filter(days.alive > threshold.day) %>%
    summarize(n.above.thres = n())
}

calculate_success_rate <- function(n.unique.antigens, n.success.antigens) {
  # Calculate how many antigens are "successful" according to a particular criteria 
  success.rate = n.success.antigens$n.above.thres/n.unique.antigens$unique.antigens
  quantile(success.rate, probs = c(0.025, .5, 0.975))*100
}

antigen_success_summary <- function(threshold.freq.levels, threshold.day.levels, antigen.frequencies) {
  # Calculate summary metrics for how a success criteria affects number of strains inlcuded 
  n.unique.antigens = count.n.antigens(antigen.frequencies)
  
  threshold.combinations = expand.grid(threshold.freq.levels, threshold.day.levels)
  colnames(threshold.combinations) = c("freq", "day")
  
  antigen.success.quantile = adply(.data = threshold.combinations, .margins = 1, function(thres) {
    thres.freq = thres$freq
    thres.day = thres$day
    success.antigens = count.success.antigens(antigen.record = antigen.frequencies, threshold.freq = thres.freq, threshold.day =  thres.day)
    quantiles = calculate.success.rate(n.unique.antigens, n.success.antigens = success.antigens)
  })
 
   num.antigen.strains = adply(.data = threshold.combinations, .margins = 1, function(thres) {
    thres.freq = thres$freq
    thres.day = thres$day
    success.antigens = count.success.antigens(antigen.record = antigen.frequencies, threshold.freq = thres.freq, threshold.day =  thres.day)
    success.antigens %>%
      summarize(min = quantile(n.above.thres, .025),
                median = quantile(n.above.thres, .5),
                max = quantile(n.above.thres, .975))
  })
  antigen.success.summary.df = left_join(antigen.success.quantile, num.antigen.strains)
  #antigen.success.summary.df$X1 = threshold.levels 
  #colnames(antigen.success.summary.df)[1] = "threshold"
  return(antigen.success.summary.df)
}

calculate_life_spans_success <- function(ant.freq.success.l, thres) {
  # takes each antigen and calculates life time above a 
  # threshold that it was above 10% 
  ant.freq.success.l %>%
    group_by(.id, antigentype) %>%
    filter(frequency > thres) %>%
    summarize(days.alive = day[n()] - day[1]) %>%
    mutate(years = days.alive/365)
}

calculate_total_life_id <- function(antigen.frequencies) {
  # This one doesn't have a threshold criteria
   antigen.frequencies %>%
    group_by(.id,antigentype) %>%
    summarize(day.emerge = day[1],
              last.day = tail(day)[6]) -> birth.death.days 
   
  birth.death.days %>% 
    mutate(life.length = ifelse(is.na(last.day), 0, last.day-day.emerge)) %>%
    dplyr::select(.id,antigentype, life.length)
}





calculate_total_life <- function(antigen.frequencies) {
 antigen.frequencies %>%
    group_by(antigentype) %>%
    summarize(day.emerge = day[1],
              last.day = last(day)) -> birth.death.days 
  
  birth.death.days %>% 
    mutate(life.length = ifelse(is.na(last.day), 0, last.day-day.emerge)) %>%
    dplyr::select(antigentype, life.length)
}


calculate_days_above_thres <- function(antigen.frequencies, threshold) {
 # Calculates the number of days an antigen type is present above a certain frequency
  antigen.frequencies %>%
    group_by(.id, antigentype) %>%
    dplyr::filter(frequency > threshold) %>%
    summarize(day.emerge = day[1],
              last.day = tail(day)[1],
              frequency.emerge = frequency[1],
              frequency.down = tail(frequency)[1]) %>%
    mutate(days.above = last.day-day.emerge) %>%
    dplyr::select(.id, antigentype, days.above) %>%
    mutate_at("antigentype", as.character) 
}

calculate_max_infected <- function(timeseries){
  # calculate the 
  timeseries %>%
    group_by(.id) %>%
    summarize(max.I = max(totalI),
              min.I = min(totalI))
}

normalize_infection <- function(meta.data, infected.range) {
  meta.data %>% 
    left_join(infected.range) %>%
    mutate(normalize.I = (infected-min.I)/(max.I-min.I)) 
}


calculate_age_of_dominant <- function(data) {
  data %>%
    group_by(.id) %>%
    dplyr::select(.id, day, postAntigen,dominant.type) -> dominant.on.day.of.emergence
  
    age.of.dominant = ddply(.data = dominant.on.day.of.emergence, .variables = ".id", function(sim) {
      for(i in 1:nrow(sim)){
        
        dominant.type = sim$dominant.type[i]
        if(dominant.type == 0) {
          sim$age[i] = (as.numeric(sim$day[i])-0)/365
        }
        else {
          day.emerge = as.numeric(sim$day[which(sim$postAntigen == dominant.type)])
          sim$age[i] = (as.numeric(sim$day[i])-day.emerge)/365
        }
      }
    return(sim)
  })
  
 data %>%
   left_join(age.of.dominant) 
}

return_success_types <- function(meta.data) {
  meta.data  %>%
    group_by(.id) %>%
    filter(success == "yes") %>%
    dplyr::select(.id, postAntigen)
}

filter_frequencies_success <- function(success.types, antigen.frequencies) {
  # Filter antigen frequencies to only the ones that are succesful
  # Works for multiple trials 
  antigen.freq.success.df = ddply(.data = success.types, .variables = ".id", function(sim) {
    types = sim$postAntigen 
    types = c(0, types)
    antigen.frequencies %>%
      filter(.id == sim$.id[1]) %>%
      filter(antigentype %in% types) -> antigen.freq.sim
    return(antigen.freq.sim)
  })
  return(antigen.freq.success.df)
}
 

vif.mer <- function (fit) {
  ## adapted from rms::vif
  
  v <- vcov(fit)
  nam <- names(fixef(fit))
  
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}


########################## Model Analysis Functions
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

##### Slicing Functions
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
  browser()
  pop.level.list = map(sample.times, collect_poplevel, trial.name)
  pop.level.df = do.call("rbind", pop.level.list) %>% arrange(day)
  
  
  ### Identify which antigens are present at these times 
  antigen.present.list = map(sample.times, find_antigens_present, trial.name)
  antigen.present.df = do.call("rbind", antigen.present.list) %>% arrange(day)
  
  
  ### Filter out those that are above the threshold of capture (1%) have already been successful
  antigen.present.df %>%
    filter(frequency > .01) -> antigen.present.captured
  
  ### Filter out ones that have reached their peak already 
  if(time.type == "first") { 
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
  } else {
    antigen.present.analyze = antigen.present.captured
  }
  
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


create_trial_dataset <- function(trial.name, sample.times) { 
  pop.level.list = map(sample.times, collect_poplevel, trial.name)
  pop.level.df = do.call("rbind", pop.level.list) %>% arrange(day)
  
  
  ### Identify which antigens are present at these times 
  antigen.present.list = map(sample.times, find_antigens_present, trial.name)
  antigen.present.df = do.call("rbind", antigen.present.list) %>% arrange(day)
  
  
  ### Filter out those that are above the threshold of capture (1%) have already been successful
  antigen.present.df %>%
    filter(frequency > .01 & frequency < .2) -> antigen.present.captured
  
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


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

#find_successful_types <- function(frequencies, threshold, length) {
  # determines succcess if a 
 # frequencies %>%
#    filter(frequency > threshold) %>%
#    group_by(antigentype) %>%
#    summarize(n = n(), 
#              max.freq  = max(frequency)) %>%
#    filter(n > length) -> successful.types
#    successful.types.id = successful.types$antigentype

  #  return(successful.types.id)
#}

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
  ### Combines all the output files for novel antigens
  # reading in and cleaning up the console file 
  print(sim.dir)
 # browser()
  console.file <- read.table(paste0(sim.dir, "/out.console.txt"), header = TRUE, fill = TRUE)
  console.file = clean_console(console.file) 
  
  # find the time when novel antigens emerge
  novel.types = find_antigen_emergence(console.file)
  
  # read in output file that tracks state of simulation when antigens emerge
  track.antigen <- read.table(paste0(sim.dir, "/out.trackAntigenSeries.txt"), header = TRUE)
  track.antigen$day = as.character(track.antigen$day)
  novel.types$day = as.character(novel.types$day)
  meta.data  = left_join(x = novel.types, y = track.antigen, "day") 
  meta.data$day = as.character(meta.data$day)
  
  # record the exact days when new novel types are generated 
  days.of.emergence = meta.data$day
  
  ## Read in viral fitness, dplyr::select on the days that correspond to emergence 
  viral.fitness <- read.table(paste0(sim.dir, "/out.viralFitnessSeries.txt"), header = TRUE)

  viral.fitness %>%
    filter(day %in% days.of.emergence) %>%
    mutate(day = as.character(day)) %>%
    dplyr::select(-simDay) %>%
    filter(!duplicated(day)) -> viral.fitness.emergence
  meta.data %>%
    left_join(viral.fitness.emergence, by = "day") -> meta.data
  
  ## combine dominant frequency circulating at time of emergence
  antigen.frequencies <- read.table(paste0(sim.dir, "/out.antigenFrequencies.txt"), header = TRUE)
  dominant.types <- find_dominant_types_at_emerge(antigen.frequencies)
 
  dominant.types %>%
    filter(day %in% days.of.emergence) -> dominant.types
  colnames(dominant.types)[2] = "dominant.type"; colnames(dominant.types)[3] = "dominant.freq"
  dominant.types$day = as.character(dominant.types$day)
  meta.data %>%
    left_join(dominant.types, by = "day") -> meta.data
  
  # combine maximum frequency the strain itself ever achieved
  maximum.freq.type <- find_max_frequency(antigen.frequencies)
  maximum.freq.type$antigentype = as.factor(maximum.freq.type$antigentype)
  meta.data$postAntigen = as.factor(meta.data$postAntigen)
  meta.data %>% left_join(maximum.freq.type, by = c("postAntigen" = "antigentype")) -> meta.data
  
  life.span <- calculate_total_life(antigen.frequencies)
  life.span$antigentype = as.factor(life.span$antigentype)
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

count_n_antigens <- function(antigen.record) {
  # Return number of unique antigens in a record for each trial
  antigen.record %>%
    group_by(.id) %>%
    summarize(unique.antigens = n_distinct(antigentype)) -> n.antigens
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
    dplyr::select(.id, antigentype, days.above) -> antigen.frequencies
  antigen.frequencies$antigentype = as.character(antigen.frequencies$antigentype)
  return(antigen.frequencies)
}

calculate_max_infected <- function(timeseries){
  # calculate the 
  timeseries %>%
    group_by(.id) %>%
    summarize(max.I = max(totalI),
              min.I = min(totalI)) -> timeseries
}

normalize_infection <- function(meta.data, infected.range) {
  meta.data %>% 
    left_join(infected.range) %>%
    mutate(normalize.I = (infected-min.I)/(max.I-min.I)) -> meta.data
  return(meta.data)
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

kappa.mer <- function (fit,
                       scale = TRUE, center = FALSE,
                       add.intercept = TRUE,
                       exact = FALSE) {
  X <- fit@pp$X
  nam <- names(fixef(fit))
  
  ## exclude intercepts
  nrp <- sum(1 * (nam == "(Intercept)"))
  if (nrp > 0) {
    X <- X[, -(1:nrp), drop = FALSE]
    nam <- nam[-(1:nrp)]
  }
  
  if (add.intercept) {
    X <- cbind(rep(1), scale(X, scale = scale, center = center))
    kappa(X, exact = exact)
  } else {
    kappa(scale(X, scale = scale, center = scale), exact = exact)
  }
}

colldiag.mer <- function (fit,
                          scale = TRUE, center = FALSE,
                          add.intercept = TRUE) {
  ## adapted from perturb::colldiag, method in Belsley, Kuh, and
  ## Welsch (1980).  look for a high condition index (> 30) with
  ## more than one high variance propotion.  see ?colldiag for more
  ## tips.
  result <- NULL
  if (center) 
    add.intercept <- FALSE
  if (is.matrix(fit) || is.data.frame(fit)) {
    X <- as.matrix(fit)
    nms <- colnames(fit)
  }
  else if (class(fit) == "mer") {
    nms <- names(fixef(fit))
    X <- fit@X
    if (any(grepl("(Intercept)", nms))) {
      add.intercept <- FALSE
    }
  }
  X <- X[!is.na(apply(X, 1, all)), ]
  
  if (add.intercept) {
    X <- cbind(1, X)
    colnames(X)[1] <- "(Intercept)"
  }
  X <- scale(X, scale = scale, center = center)
  
  svdX <- svd(X)
  svdX$d
  condindx <- max(svdX$d)/svdX$d
  dim(condindx) <- c(length(condindx), 1)
  
  Phi = svdX$v %*% diag(1/svdX$d)
  Phi <- t(Phi^2)
  pi <- prop.table(Phi, 2)
  colnames(condindx) <- "cond.index"
  if (!is.null(nms)) {
    rownames(condindx) <- nms
    colnames(pi) <- nms
    rownames(pi) <- nms
  } else {
    rownames(condindx) <- 1:length(condindx)
    colnames(pi) <- 1:ncol(pi)
    rownames(pi) <- 1:nrow(pi)
  }         
  
  result <- data.frame(cbind(condindx, pi))
  zapsmall(result)
}

maxcorr.mer <- function (fit,
                         exclude.intercept = TRUE) {
  so <- summary(fit)
  corF <- so@vcov@factors$correlation
  nam <- names(fixef(fit))
  
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0 & exclude.intercept) {
    corF <- corF[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  corF[!lower.tri(corF)] <- 0
  maxCor <- max(corF)
  minCor <- min(corF)
  if (abs(maxCor) > abs(minCor)) {
    zapsmall(maxCor)
  } else {
    zapsmall(minCor)
  }
}